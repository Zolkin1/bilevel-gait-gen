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
                                 mpc_(info, robot_urdf) {

        // TODO: Get this from IC
        std::array<std::array<double, 3>, 4> ee_pos{};
        ee_pos.at(0) = {0.2, 0.2, 0};
        ee_pos.at(1) = {0.2, -0.2, 0};
        ee_pos.at(2) = {-0.2, 0.2, 0};
        ee_pos.at(3) = {-0.2, -0.2, 0};

        mpc_.SetStateTrajectoryWarmStart(warm_start_states);

        mpc_.SetDefaultGaitTrajectory(mpc::Gaits::Trot, 2, ee_pos);

        matrix_t Q = matrix_t::Zero(24, 24);
        Q.topLeftCorner<6,6>() = matrix_t::Zero(6,6);
        Q(6,6) = 30;
        Q(7,7) = 30;
        Q(8,8) = 30;

        const vector_t des_alg = mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(state_des, warm_start_states.at(0));

        std::cout << "Desired state: \n" << des_alg << std::endl;

        mpc_.AddQuadraticTrackingCost(des_alg, Q);
        mpc_.AddForceCost(0.01);
        mpc_.SetQuadraticFinalCost(50*Q);
        mpc_.SetLinearFinalCost(-50*Q*des_alg);
        prev_time_ = 0;
        computed_ = false;
    }

    vector_t MPCController::ComputeControlAction(const vector_t& q, const vector_t& v,
                                                 const vector_t& a, const controller::Contact& contact,
                                                 double time) { //, double time) {
        vector_t state = ReconstructState(q, v, a);
//        if ((time - prev_time_ > 0.1 || time == 0)) {
        if (!computed_) {
            for (int i = 0; i < 1; i++) {
                mpc_.Solve(state, time);
            }
            prev_time_ = time;
            computed_ = true;
            mpc_.PrintStats();
        }
        if (time < 1) {
            qp_controller_.UpdateTargetConfig(mpc_.GetTargetConfig(time));
            qp_controller_.UpdateTargetVel(mpc_.GetTargetVelocity(time));
            qp_controller_.UpdateTargetAcc(mpc_.GetTargetAcc(time));
//          qp_controller_.UpdateTargetForce(mpc_.GetTargetForce());
            qp_controller_.UpdateDesiredContacts(mpc_.GetDesiredContacts(time));
        }

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

} // controller