//
// Created by zolkin on 11/29/23.
//

#include "mpc_controller.h"

namespace controller {
    MPCController::MPCController(double control_rate, std::string robot_urdf, const std::string &foot_type, int nv,
                                 const Eigen::VectorXd &torque_bounds, double friction_coef,
                                 const std::vector<double>& base_pos_gains, const std::vector<double>& base_ang_gains,
                                 const std::vector<double>& joint_gains, double leg_weight, double torso_weight,
                                 double force_weight, mpc::MPCInfo info,
                                 const std::vector<vector_t>& warm_start_states,
                                 const vector_t& state_des) : Controller(control_rate, robot_urdf, foot_type),
                                 qp_controller_(control_rate, robot_urdf, foot_type, nv,
                                               torque_bounds, friction_coef,
                                               base_pos_gains,
                                               base_ang_gains,
                                               joint_gains,
                                               leg_weight,
                                               torso_weight,
                                               force_weight,
                                               4),
                                 mpc_(info, robot_urdf),
                                  traj_(10, 24, 12, mpc_.CreateDefaultSwitchingTimes(2,4,1.0),
                                        0.015, 0.75, 0.0) {

        // TODO: Get this from IC
        std::array<std::array<double, 3>, 4> ee_pos{};
        ee_pos.at(0) = {0.2, 0.2, 0};
        ee_pos.at(1) = {0.2, -0.2, 0};
        ee_pos.at(2) = {-0.2, 0.2, 0};
        ee_pos.at(3) = {-0.2, -0.2, 0};

        mpc_.SetStateTrajectoryWarmStart(warm_start_states);

        mpc_.SetDefaultGaitTrajectory(mpc::Gaits::Trot, 3, ee_pos); // TODO: Get this from the config file

        matrix_t Q = 2*matrix_t::Identity(24, 24);
        Q.topLeftCorner<6,6>() = matrix_t::Zero(6, 6);
        Q(0,0) = 0.0;
        Q(1,1) = 0.0;
        Q(2,2) = 10;
        Q(6,6) = 300; //300;
        Q(7,7) = 300; //300;
        Q(8,8) = 750; //450;
        Q(9,9) = 200;
        Q(10,10) = 200;
        Q(11,11) = 200;

        Q(12,12) = 600;
        Q(15,15) = 600;
        Q(18,18) = 600;
        Q(21,21) = 600;

        Q(13,13) = 50;
        Q(16,16) = 50;
        Q(19,19) = 50;
        Q(22,22) = 50;

        Q(14,14) = 10;
        Q(17,17) = 10;
        Q(20,20) = 10;
        Q(23,23) = 10;

        const vector_t des_alg = mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(state_des, warm_start_states.at(0));

        std::cout << "Desired state: \n" << des_alg << std::endl;

        mpc_.AddQuadraticTrackingCost(des_alg, Q);
        mpc_.AddForceCost(0.00);
        mpc_.SetQuadraticFinalCost(500*Q);
        mpc_.SetLinearFinalCost(-500*Q*des_alg);

        prev_time_ = 0;
        computed_ = false;
    }

    void MPCController::InitSolver(const vector_t& state) {
        mpc_.CreateInitialRun(state);
        mpc_.PrintStats();

        state_ = state;
        time_ = 0;
        q_des_ = mpc_.GetTargetConfig(time_);
        v_des_ = mpc_.GetTargetVelocity(time_);
        a_des_ = mpc_.GetTargetAcc(time_);
        force_des_ = mpc_.GetForceTarget(time_);
        traj_ = mpc_.GetTrajectory();
        mpc_computations_ = std::thread(&MPCController::MPCUpdate, this);
    }

    vector_t MPCController::ComputeControlAction(const vector_t& q, const vector_t& v,
                                                 const vector_t& a, const controller::Contact& contact,
                                                 double time) {
        // TODO: Investigate all the conversions - may be causing an issue
        const vector_t state = ReconstructState(q, v, a);
        state_time_mut_.lock();
        state_ = state;
        time_ = time;
        state_time_mut_.unlock();

//        if (time - prev_time_ >= 0.015) {
//            mpc_.GetRealTimeUpdate(50, state, time);
//            prev_time_ = time;
//            std::cout << "movement error at " << time << "s \n" << traj_.GetState(1) - state << std::endl;
//        }
//        traj_ = mpc_.GetTrajectory();

        mpc_res_mut_.lock();
        int node = mpc_.GetNode(time)+1;

        qp_controller_.UpdateTargetConfig(traj_.GetState(node).tail(num_inputs_ + FLOATING_BASE_OFFSET));
//        qp_controller_.UpdateTargetVel(traj_.GetInputs().GetVels(time).tail(num_inputs_ + FLOATING_BASE_OFFSET)); // TODO: Fix
//        qp_controller_.UpdateTargetAcc(a_des_);
        qp_controller_.UpdateForceTargets(force_des_);
        mpc_res_mut_.unlock();

        Contact contact1 = mpc_.GetDesiredContacts(time);
//        contact1.in_contact_ = {false, true, true, false};
//        contact1.contact_frames_ = {14, 24, 34, 44};

        qp_controller_.UpdateDesiredContacts(contact1);

        return qp_controller_.ComputeControlAction(q, v, a, contact, time);
    }

    vector_t MPCController::ReconstructState(const controller::vector_t& q, const controller::vector_t& v,
                                             const controller::vector_t& a) const {

        vector_t state = vector_t::Zero(q.size() + 6);
        state.head<3>() = v.head<3>() * mpc_.GetModel().GetMass();

        state.segment<3>(3) = v.segment<3>(3) * mpc_.GetModel().GetMass();

        Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(q.segment<4>(3)));
        // Note the warning on the pinocchio function!
        pinocchio::quaternion::firstOrderNormalize(quat);
        state(9) = quat.x();
        state(10) = quat.y();
        state(11) = quat.z();
        state(12) = quat.w();

        state.segment<mpc::CentroidalModel::POS_VARS>(mpc::CentroidalModel::MOMENTUM_OFFSET) = q.head<mpc::CentroidalModel::POS_VARS>();

        state.tail(q.size() - mpc::CentroidalModel::FLOATING_VEL_OFFSET) =
                q.tail(q.size() - mpc::CentroidalModel::FLOATING_VEL_OFFSET);

        return state;
    }

    void MPCController::MPCUpdate() {
        vector_t state;
        double time;

        while(true) {   // TODO: Add an end condition
            state_time_mut_.lock();
            state = state_;        // TODO: Make sure to get a copy, not a reference
            time = time_;
            state_time_mut_.unlock();

            mpc_.GetRealTimeUpdate(200, state, time);

            //            std::this_thread::sleep_for(std::chrono::milliseconds(40));

            mpc_res_mut_.lock();
//            q_des_ = mpc_.GetNextTargetConfig();
//            v_des_ = mpc_.GetTargetVelocity(time);
//            a_des_ = mpc_.GetTargetAcc(time);
//            force_des_ = mpc_.GetForceTarget(time);
            traj_ = mpc_.GetTrajectory();
            mpc_res_mut_.unlock();
        }
    }

} // controller