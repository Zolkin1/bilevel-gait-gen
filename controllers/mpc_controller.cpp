//
// Created by zolkin on 11/29/23.
//

#include "mpc_controller.h"
#include "models/single_rigid_body_model.h"

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/frames.hpp"

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

        num_polys_ = num_polys;

        // TODO: Do I need the model if the basic controller also has pinocchio data?

        force_des_.resize(4);

        mpc_.SetStateTrajectoryWarmStart(warm_start_states);


        const vector_t des_alg = model_.ConvertManifoldStateToTangentState(state_des, warm_start_states.at(0));

        std::cout << "Desired state: \n" << des_alg << std::endl;

        mpc_.AddQuadraticTrackingCost(des_alg, Q);
        mpc_.AddForceCost(0.00);
        mpc_.SetQuadraticFinalCost(1*Q);
        mpc_.SetLinearFinalCost(-1*Q*des_alg);

        kp_joints_ = kp_joint_gains;
        kv_joints_ = kd_joint_gains;

        prev_time_ = 0;
        computed_ = false;

        run_num = 0;

        for (auto& ee : fk_traj_) {
            ee.resize(info.num_nodes + 1);
        }
    }

    void MPCController::InitSolver(const vector_t& full_body_state, const vector_t& mpc_state) {
        state_ = mpc_state;
        full_body_state_ = full_body_state;
        q_des_ = full_body_state;
        ee_locations_ = model_.GetEndEffectorLocations(full_body_state);
        mpc_.SetDefaultGaitTrajectory(mpc::Gaits::Trot, num_polys_, ee_locations_);

//        std::array<std::array<double, 3>, 4> ee_pos{};
//        ee_pos.at(0) = {0.2, 0.2, 0};
//        ee_pos.at(1) = {0.2, -0.2, 0};
//        ee_pos.at(2) = {-0.2, 0.2, 0};
//        ee_pos.at(3) = {-0.2, -0.2, 0};
//
//        ee_locations_.resize(4);
//        for (int i = 0; i < ee_locations_.size(); i++) {
//            for (int j = 0; j < 3; j++) {
//                ee_locations_.at(i)(j) = ee_pos.at(i).at(j);
//            }
//        }

        mpc_.CreateInitialRun(mpc_state, ee_locations_);
        mpc_.PrintStats();

//        traj_viz_mut_.lock();
        traj_ = mpc_.GetTrajectory();
        FullBodyTrajUpdate(traj_);
//        traj_viz_mut_.unlock();

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

        if (time > time_) {
            time_ = time;
        }

//        if (time > 1) {
//            time = 1;
//        }

        ee_locations_ = ee_locations;
        state_time_mut_.unlock();

        // For now, single thread everything

        // Get MPC update (not gait opt yet)
        mpc::Trajectory traj = traj_;
//        if (!(run_num % 10)) {
//            mpc_.AdjustForCurrentContacts(time, contact);
//            traj = mpc_.GetRealTimeUpdate(state, time, ee_locations, false);
//        }
        FullBodyTrajUpdate(traj);
        UpdateTrajViz();

//        mpc_res_mut_.lock();
        GetTargetsFromTraj(traj, time); // TODO: IK for velocties.
        // TODO: Check GRF
        Contact contact1 = traj.GetDesiredContacts(time);  // TODO: Implement
        contact1.contact_frames_ = contact.contact_frames_;
//        mpc_res_mut_.unlock();

//        q_des_ << 0., 0., 0.36, 0, 0, 0, 1.0, -0.22, -0.2, -0.977,
//                0.02, 0.2, -0.977,  0.02, 0.45, -1.03, -0.02, 0.45, -1.03;
//        v_des_.setZero();

        vector_t tau_fb = kp_joints_.cwiseProduct((q_des_ - q).tail(num_inputs_))
                + kv_joints_.cwiseProduct((v_des_ - v).tail(num_inputs_));
//        tau_fb.setZero();

        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);
        // Get C
        pinocchio::computeCoriolisMatrix(pin_model_, *pin_data_, q, v);
        // Get g
        pinocchio::computeGeneralizedGravity(pin_model_, *pin_data_, q);
        // Get M
        pinocchio::crba(pin_model_, *pin_data_,  q);
        // Make M symmetric
        pin_data_->M.triangularView<Eigen::StrictlyLower>() =
                pin_data_->M.transpose().triangularView<Eigen::StrictlyLower>();

        const int num_sw = (4-contact.GetNumContacts())*3;
        const int num_st = contact.GetNumContacts()*3;
        vector_t tau_ff(num_inputs_); //pin_data_->g.tail(num_inputs_);

        for (int ee = 0; ee < contact1.in_contact_.size(); ee++) {
            if (contact1.in_contact_.at(ee)) {
                // Stance foot
//                const mpc::matrix_33t J = model_.GetEEJacobian(ee, q);
//                const mpc::matrix_33t Jdot = model_.GetEEJacobianDeriv(ee, q, v);

                matrix_t J = matrix_t::Zero(6, 18);
                pinocchio::computeFrameJacobian(pin_model_, *pin_data_, q, contact1.contact_frames_.at(ee),
                                                pinocchio::LOCAL_WORLD_ALIGNED, J); // TODO: Consider frame

                // force in generalized coordinates
                vector_t qf = -J.transpose().leftCols<3>()*traj.GetForce(ee, time);

//                std::cout << "J^T: \n" << J.transpose() << std::endl;
//                std::cout << "force: \n" << traj_.GetForce(ee, time) << std::endl;
//
//                std::cout << "qf: \n" << qf << std::endl;

                tau_ff.segment(3*ee, 3) = qf.segment<3>(6 + 3*ee); // .setZero();
            } else {
                // Swing foot
                tau_ff.segment(3*ee, 3) = (-pin_data_->g + pin_data_->C*v + pin_data_->M*a).segment(6 + 3*ee,3);
            }
//            std::cout << "ee: " << ee << ", contact: " << contact.in_contact_.at(ee) << std::endl;
        }
//        tau_ff.setZero();

        vector_t control_action = vector_t::Zero(3*num_inputs_);    // joints, joint velocities, torques

        control_action.head(num_inputs_) = q_des_.tail(num_inputs_);
        control_action.segment(num_inputs_, num_inputs_) = vector_t::Zero(12); //v_des_.tail(num_inputs_);
        control_action.tail(num_inputs_) = tau_fb + tau_ff;

//        for (int ee = 0; ee < contact1.in_contact_.size(); ee++) {
            // TODO: Should I be going on actual or desired contacts?
//            if (contact1.in_contact_.at(ee)) {
//                control_action.segment(num_inputs_*2 + 3*ee, 3) = ComputeGroundForceTorques(ee); // TODO: Make not hard coded
//            } else {
//                control_action.segment(num_inputs_*2 + 3*ee, 3) = ComputeSwingTorques(ee); // TODO: Make not hard coded
//            }
//        }

//        control_action.segment(num_inputs_*2 + 6, 3) = vector_t::Zero(3);

//        return control_action;

        qp_controller_.UpdateTargetConfig(q_des_);
        qp_controller_.UpdateTargetVel(v_des_);
        vector_t force_des(contact1.GetNumContacts()*3);
        int j = 0;
        for (int i = 0; i < contact1.in_contact_.size(); i++) {
            if (contact1.in_contact_.at(i)) {
                force_des.segment<3>(3*j) = traj.GetForce(i, time);
                j++;
            }
        }
        qp_controller_.UpdateForceTargets(force_des);

        qp_controller_.UpdateDesiredContacts(contact1);

        traj_ = traj;

        run_num++;

        return qp_controller_.ComputeControlAction(q, v, a, contact1, time);
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
        Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(q.segment<4>(3)));
        // Note the warning on the pinocchio function!
        pinocchio::quaternion::firstOrderNormalize(quat);
        state(mpc::SingleRigidBodyModel::ORIENTATION_START) = quat.x();
        state(mpc::SingleRigidBodyModel::ORIENTATION_START + 1) = quat.y();
        state(mpc::SingleRigidBodyModel::ORIENTATION_START + 2) = quat.z();
        state(mpc::SingleRigidBodyModel::ORIENTATION_START + 3) = quat.w();

        // Angular Momentum
        // TODO: Is the angular velocity correct?
        // TODO: The state is the angular momentum, not just the velocity
        state.segment<POS_VARS>(mpc::SingleRigidBodyModel::ANG_VEL_START + 1) =
                model_.GetIr() * v.segment<POS_VARS>(mpc::SingleRigidBodyModel::QUAT_SIZE);

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
//        vector_t state;
//        double time = -1;
//        std::vector<mpc::vector_3t> ee_locations(4);
//
//        while(true) {
//            state_time_mut_.lock();
//            double time_act = time_;
//            state_time_mut_.unlock();
//            if (time_act - time >= 0.1) {
//                state_time_mut_.lock();
//                state = state_;
//                time = time_;
//                ee_locations = ee_locations_;
//                state_time_mut_.unlock();
//
//                mpc::Trajectory traj = mpc_.GetRealTimeUpdate(state, time, ee_locations, false);
//
//                // ----------------------- IK On Trajectory ----------------------- //
//                FullBodyTrajUpdate(traj);
//                // ----------------------- ----------------------- //
//
//                mpc_res_mut_.lock();
//                traj_ = traj;
//                UpdateTrajViz();
//                mpc_res_mut_.unlock();

//                gait_optimizer_.UpdateSizes(mpc_.GetNumDecisionVars(), mpc_.GetNumConstraints());
//                mpc_.ComputeDerivativeTerms();
//                mpc_.GetQPPartials(gait_optimizer_.GetQPPartials());

//                gait_optimizer_.SetParameterPartials()

//                gait_optimizer_.ComputeCostFcnDerivWrtContactTimes();

//                gait_optimizer_.OptimizeContactTimes();

//                mpc_.UpdateContactTimes(gait_optimizer_.GetContactTimes());

//            }
//        }
    }

    std::vector<std::vector<Eigen::Vector3d>> MPCController::GetTrajViz() {
//        traj_viz_mut_.lock();
        std::vector<std::vector<Eigen::Vector3d>> fk_traj = fk_traj_;
//        traj_viz_mut_.unlock();
        return fk_traj;
    }

    void MPCController::UpdateTrajViz() {
//        traj_viz_mut_.lock();
        fk_traj_ = mpc_.CreateVizData();
//        traj_viz_mut_.unlock();
    }

    void MPCController::FullBodyTrajUpdate(mpc::Trajectory& traj) {
        vector_t state_guess = q_des_;
        for (int i = 0; i < traj.GetStates().size(); i++) {
            std::vector<mpc::vector_3t> traj_ee_locations(4);
            for (int ee = 0; ee < 4; ee++) {
                // Grab state and end effector locations at that time
                traj_ee_locations.at(ee) = traj.GetEndEffectorLocation(ee, traj.GetTime(i));
            }
            traj.UpdateFullConfig(i, model_.InverseKinematics(traj.GetState(i),
                                                              traj_ee_locations, state_guess,
                                                              info_.joint_bounds_ub,
                                                              info_.joint_bounds_lb));
            state_guess = traj.GetFullConfig(i);

            if (i >= 1) {
                // TODO: Do the velocity IK (better vel estimate...)

                vector_t vel(model_.GetFullModelConfigSpace() - 1);
                vel.head<6>() << traj.GetState(i).segment<3>(3), traj.GetState(i).segment<3>(9);
                vel.tail(num_inputs_) = (traj.GetFullConfig(i) - traj.GetFullConfig(i-1)).tail(num_inputs_)/info_.integrator_dt;

                traj.UpdateFullVelocity(i-1, vel);
            }
        }
        traj.UpdateFullVelocity(traj.GetStates().size() - 1, vector_t::Zero(pin_model_.nv));
    }

    void MPCController::GetTargetsFromTraj(const mpc::Trajectory& traj, double time) {
        if (time < traj.GetTime(0)) {
            time = traj.GetTime(0);
        }

        int node = traj.GetNode(time)+1; // TODO: +1 or +0

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

        torques = J.transpose()*force_des_.at(end_effector); // TODO: Investigate R

        return torques;
    }

    vector_t MPCController::ComputeSwingTorques(int end_effector) {
        mpc::vector_3t acc;
        acc.setZero();

        for (int ee = 0; ee < 4; ee++) {
            acc += force_des_.at(ee);
        }
        acc = acc/model_.GetMass();


        const mpc::matrix_33t Lambda = model_.GetOperationalSpaceInertia(end_effector, q_des_);
        const vector_t NL = model_.GetNonlinearEffects(q_des_, v_des_);

        const mpc::matrix_33t J = model_.GetEEJacobian(end_effector, q_des_);
        const mpc::matrix_33t Jdot = model_.GetEEJacobianDeriv(end_effector, q_des_, v_des_);

        vector_t torques = J.transpose()*Lambda*(acc - Jdot*v_des_.segment<3>(6 + 3*end_effector))
                + NL.segment<3>(6 + 3*end_effector);
                //(C*v_des_).segment<3>(6 + 3*end_effector) + G.segment<3>(6 + 3*end_effector);   // TODO: Check the indexing

        return torques;
    }

    std::vector<Eigen::Vector2d> MPCController::GetEEBoxCenter() {
        // TODO: Deal with multithreading
        return mpc_.GetEEBoxCenter();
    }

} // controller