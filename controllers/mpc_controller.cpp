//
// Created by zolkin on 11/29/23.
//

#include "mpc_controller.h"
#include "single_rigid_body_model.h"

namespace controller {
    MPCController::MPCController(double control_rate, std::string robot_urdf, const std::string &foot_type, int nv,
                                 const Eigen::VectorXd &torque_bounds, double friction_coef,
                                 const std::vector<double>& base_pos_gains, const std::vector<double>& base_ang_gains,
                                 const vector_t& kp_joint_gains,
                                 const vector_t& kd_joint_gains,
                                 double leg_weight, double torso_weight,
                                 double force_weight, mpc::MPCInfo info,
                                 const std::vector<vector_t>& warm_start_states,
                                 const vector_t& state_des,
                                 int num_polys,
                                 const matrix_t& Q) : Controller(control_rate, robot_urdf, foot_type),
                                 qp_controller_(control_rate, robot_urdf, foot_type, nv,
                                               torque_bounds, friction_coef,
                                               base_pos_gains,
                                               base_ang_gains,
                                               kp_joint_gains,
                                               kd_joint_gains,
                                               leg_weight,
                                               torso_weight,
                                               force_weight,
                                               4),
                                 mpc_(info, robot_urdf),
                                 traj_(10, 24, 12, mpc_.CreateDefaultSwitchingTimes(2,4,1.0),
                                        0.015, 0.75, 0.0),
                                 fk_traj_(5),
                                 gait_optimizer_(4, 10, 10, 10, 1, 0.05),
                                 info_(info),
                                 model_(robot_urdf, info.ee_frames, info.discretization_steps, info.integrator_dt, info.nom_state){  // TODO: Not hard coded

        // TODO: Do I need the model if the basic controller also has pinocchio data?

        force_des_.resize(4);

        // TODO: Get this from IC
        std::array<std::array<double, 3>, 4> ee_pos{};
        ee_pos.at(0) = {0.2, 0.2, 0};
        ee_pos.at(1) = {0.2, -0.2, 0};
        ee_pos.at(2) = {-0.2, 0.2, 0};
        ee_pos.at(3) = {-0.2, -0.2, 0};

        mpc_.SetStateTrajectoryWarmStart(warm_start_states);

        mpc_.SetDefaultGaitTrajectory(mpc::Gaits::Trot, num_polys, ee_pos);

        const vector_t des_alg = model_.ConvertManifoldStateToTangentState(state_des, warm_start_states.at(0));

        std::cout << "Desired state: \n" << des_alg << std::endl;

        mpc_.AddQuadraticTrackingCost(des_alg, Q);
        mpc_.AddForceCost(0.00);
        mpc_.SetQuadraticFinalCost(1*Q);
        mpc_.SetLinearFinalCost(-1*Q*des_alg);

        prev_time_ = 0;
        computed_ = false;

        for (auto& ee : fk_traj_) {
            ee.resize(info.num_nodes + 1);
        }
    }

    void MPCController::InitSolver(const vector_t& full_body_state, const vector_t& mpc_state) {
        state_ = mpc_state;
//        ee_locations_ = model_.GetEndEffectorLocations(full_body_state);
        std::array<std::array<double, 3>, 4> ee_pos{};
        ee_pos.at(0) = {0.2, 0.2, 0};
        ee_pos.at(1) = {0.2, -0.2, 0};
        ee_pos.at(2) = {-0.2, 0.2, 0};
        ee_pos.at(3) = {-0.2, -0.2, 0};

        ee_locations_.resize(4);
        for (int i = 0; i < ee_locations_.size(); i++) {
            for (int j = 0; j < 3; j++) {
                ee_locations_.at(i)(j) = ee_pos.at(i).at(j);
            }
        }

        // TODO: Debug in release mode
        mpc_.CreateInitialRun(mpc_state, ee_locations_);
        mpc_.PrintStats();

        traj_viz_mut_.lock();
        traj_ = mpc_.GetTrajectory();
        FullBodyTrajUpdate(traj_);
        traj_viz_mut_.unlock();

        time_ = 0;

        GetTargetsFromTraj(traj_, time_);

        UpdateTrajViz();
        mpc_computations_ = std::thread(&MPCController::MPCUpdate, this);
    }

    vector_t MPCController::ComputeControlAction(const vector_t& q, const vector_t& v,
                                                 const vector_t& a, const controller::Contact& contact,
                                                 double time) {
        const vector_t state = ReconstructState(q, v, a);
        std::vector<mpc::vector_3t> ee_locations = model_.GetEndEffectorLocations(q);

        state_time_mut_.lock();
        state_ = state;
        time_ = time;
        ee_locations_ = ee_locations;
        state_time_mut_.unlock();

        mpc_res_mut_.lock();
        GetTargetsFromTraj(traj_, time_);
        mpc_res_mut_.unlock();

        // TODO: DMA
        vector_t control_action = vector_t::Zero(3*num_inputs_);    // joints, joint velocities, torques
        control_action.head(num_inputs_) = q_des_.tail(num_inputs_);
        control_action.segment(num_inputs_, num_inputs_) = v_des_.tail(num_inputs_);
        for (int ee = 0; ee < contact.in_contact_.size(); ee++) {
            if (contact.in_contact_.at(ee)) {
                control_action.segment(num_inputs_*2 + 3*ee, 3) = ComputeGroundForceTorques(ee); // TODO: Make not hard coded
            } else {
                control_action.segment(num_inputs_*2 + 3*ee, 3) = ComputeSwingTorques(ee); // TODO: Make not hard coded
            }
        }

        return control_action;

        // Recover torques using ID

//        Contact contact1 = mpc_.GetDesiredContacts(time);
//        qp_controller_.UpdateDesiredContacts(contact1);

//        return qp_controller_.ComputeControlAction(q, v, a, contact, time);
    }

    vector_t MPCController::ReconstructState(const controller::vector_t& q, const controller::vector_t& v,
                                             const controller::vector_t& a) const {

        vector_t state = vector_t::Zero(mpc_.GetModel()->GetNumManifoldStates());

        // TODO: Is q head COM or someother frame?
        // COM position
        state.head<POS_VARS>() = q.head<3>();

        // COM momentum
        state.segment<POS_VARS>(mpc::SingleRigidBodyModel::LIN_MOM_START) = v.head<3>() * mpc_.GetModel()->GetMass();

        // Orientation
        state.segment<mpc::SingleRigidBodyModel::QUAT_SIZE>(mpc::SingleRigidBodyModel::ORIENTATION_START) =
                q.segment<mpc::SingleRigidBodyModel::QUAT_SIZE>(POS_VARS);

        // Angular Momentum
        state.segment<POS_VARS>(mpc::SingleRigidBodyModel::ANG_VEL_START) =
                v.segment<POS_VARS>(mpc::SingleRigidBodyModel::QUAT_SIZE);

//        state.segment<3>(3) = v.segment<3>(3) * mpc_.GetModel()->GetMass();

//        Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(q.segment<4>(3)));
        // Note the warning on the pinocchio function!
//        pinocchio::quaternion::firstOrderNormalize(quat);
//        state(9) = quat.x();
//        state(10) = quat.y();
//        state(11) = quat.z();
//        state(12) = quat.w();
//
//        state.segment<mpc::CentroidalModel::POS_VARS>(mpc::CentroidalModel::MOMENTUM_OFFSET) = q.head<mpc::CentroidalModel::POS_VARS>();
//
//        state.tail(q.size() - mpc::CentroidalModel::FLOATING_VEL_OFFSET) =
//                q.tail(q.size() - mpc::CentroidalModel::FLOATING_VEL_OFFSET);

        return state;
    }

    vector_t MPCController::NormalizeQuat(const vector_t& state) {
        vector_t state_norm = state;
        Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(state.segment<4>(9)));
        // Note the warning on the pinocchio function!
        pinocchio::quaternion::firstOrderNormalize(quat);
        state_norm(9) = quat.x();
        state_norm(10) = quat.y();
        state_norm(11) = quat.z();
        state_norm(12) = quat.w();

        return state_norm;
    }

    void MPCController::MPCUpdate() {
        vector_t state;
        double time = -1;
        std::vector<mpc::vector_3t> ee_locations(4);

        while(true) {
            state_time_mut_.lock();
            double time_act = time_;
            state_time_mut_.unlock();
            if (time_act - time >= 0.1) {
                state_time_mut_.lock();
                state = state_;
                time = time_;
                ee_locations = ee_locations_;
                state_time_mut_.unlock();

                mpc::Trajectory traj = mpc_.GetRealTimeUpdate(state, time, ee_locations);

                // ----------------------- IK On Trajectory ----------------------- //
                FullBodyTrajUpdate(traj);
                // ----------------------- ----------------------- //

                mpc_res_mut_.lock();
                traj_ = traj;
                UpdateTrajViz();
                mpc_res_mut_.unlock();

//                gait_optimizer_.UpdateSizes(mpc_.GetNumDecisionVars(), mpc_.GetNumConstraints());
//                mpc_.ComputeDerivativeTerms();
//                mpc_.GetQPPartials(gait_optimizer_.GetPartials());

//                gait_optimizer_.SetParameterPartials()

//                gait_optimizer_.ComputeCostFcnDerivWrtContactTimes();

//                gait_optimizer_.OptimizeContactTimes();

//                mpc_.UpdateContactTimes(gait_optimizer_.GetContactTimes());

            }
        }
    }

    std::vector<std::vector<Eigen::Vector3d>> MPCController::GetTrajViz() {
        traj_viz_mut_.lock();
        std::vector<std::vector<Eigen::Vector3d>> fk_traj = fk_traj_;
        traj_viz_mut_.unlock();
        return fk_traj;
    }

    void MPCController::UpdateTrajViz() {
        traj_viz_mut_.lock();
        fk_traj_ = mpc_.CreateVizData();
        traj_viz_mut_.unlock();
    }

    void MPCController::FullBodyTrajUpdate(mpc::Trajectory& traj) {
        // TODO: Figure out a better initial guess
        vector_t prev_state = vector_t::Zero(model_.GetNumManifoldStates());
        for (int i = 0; i < 20; i++) {
            std::vector<mpc::vector_3t> traj_ee_locations(4);
            for (int ee = 0; ee < 4; ee++) {
                // Grab state and end effector locations at that time
                traj_ee_locations.at(ee) = traj.GetEndEffectorLocation(ee, traj.GetTime(i));
            }
            traj.UpdateFullConfig(i, model_.InverseKinematics(traj.GetState(i),
                                                              traj_ee_locations, prev_state,
                                                              info_.joint_bounds_ub,
                                                              info_.joint_bounds_lb));

            // TODO: Do the velocity IK
            traj.UpdateFullVelocity(i, vector_t::Zero(model_.GetFullModelConfigSpace()-1));
        }
    }

    void MPCController::GetTargetsFromTraj(const mpc::Trajectory& traj, double time) {
        const int node = traj.GetNode(time)+1;

        q_des_ = traj.GetFullConfig(node);
        v_des_ = traj.GetFullVelocity(node);

        // TODO: Do I want to use the ZoH or the full spline force?
        for (int ee = 0; ee < 4; ee++) {
            force_des_.at(ee) = traj.GetForce(ee, time + info_.integrator_dt);
        }
    }

    vector_t MPCController::ComputeGroundForceTorques(int end_effector) {
        vector_t torques(3);    // TODO: Make not hard coded

        mpc::matrix_33t J = model_.GetEEJacobian(end_effector, q_des_);
        mpc::matrix_33t R = model_.GetBaseRotationMatrix(q_des_);

        torques = J*R*force_des_.at(end_effector);

        return torques;
    }

    vector_t MPCController::ComputeSwingTorques(int end_effector) {
        mpc::vector_3t acc;
        acc.setZero();

        for (int ee = 0; ee < 4; ee++) {
            acc += force_des_.at(ee);
        }
        acc = acc/model_.GetMass();


        const mpc::matrix_33t J = model_.GetEEJacobian(end_effector, q_des_);
        const mpc::matrix_33t Jdot = model_.GetEEJacobianDeriv(end_effector, q_des_, v_des_);
        const mpc::matrix_33t Lambda = model_.GetOperationalSpaceInertia(end_effector, q_des_);

        const mpc::matrix_3t C = model_.GetCoriolisMat(q_des_, v_des_);
        const mpc::vector_3t G = model_.GetGravityVec(q_des_);

        vector_t torques = J.transpose()*Lambda*(acc - Jdot*v_des_) + C*v_des_ + G;

        return torques;
    }

} // controller