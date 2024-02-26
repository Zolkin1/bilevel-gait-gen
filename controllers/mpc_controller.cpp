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
#include "pinocchio/algorithm/center-of-mass.hpp"

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
                                 const matrix_t& Q,
                                 int gait_opt_freq) : Controller(control_rate, robot_urdf, foot_type),
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
                                 model_(robot_urdf, info.ee_frames, info.discretization_steps, info.integrator_dt, info.nom_state),
                                 gait_opt_(4, 10, 10, 10, 1, 0.05) {  // TODO: Not hard coded

        num_polys_ = num_polys;

        gait_opt_freq_ = gait_opt_freq;

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

        log_file_.open("mpc_controller_log.txt");
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

        utils::Timer timer("mpc controller");
        timer.StartTimer();

        if (kp_joints_.size() != num_inputs_) {
            throw std::runtime_error("kp joints wrong size!");
        }

        if (kv_joints_.size() != num_inputs_) {
            throw std::runtime_error("kv joints wrong size!");
        }

        const vector_t state = ReconstructState(q, v, a);
//        std::cout << "expected lin momentum: " << traj_.GetState(traj_.GetNode(time)).segment<3>(3).transpose() << std::endl;
//        std::cout << "expected ang momentum: " << traj_.GetState(traj_.GetNode(time)).tail(3).transpose() << std::endl;
        std::vector<mpc::vector_3t> ee_locations = model_.GetEndEffectorLocations(q);

        state_time_mut_.lock(); // -------- Setting MPC Variables -------- //
        state_ = state;

        if (time > time_) {
            time_ = time;
        }

        ee_locations_ = ee_locations;
        contact_ = contact;
        state_time_mut_.unlock(); // -------- End Setting MPC Variables -------- //

        // For now, single thread everything
        // Get MPC update (not gait opt yet)
//        mpc::Trajectory traj = traj_;
//        if (!(run_num % 10)) {
//            mpc_.AdjustForCurrentContacts(time, contact);
//            traj = mpc_.GetRealTimeUpdate(state, time, ee_locations, false);
//        }
//        FullBodyTrajUpdate(traj);
//        UpdateTrajViz();

        mpc_res_mut_.lock(); // -------- Accessing MPC Result -------- //
        GetTargetsFromTraj(traj_, time);
        Contact contact1 = traj_.GetDesiredContacts(time);
        contact1.contact_frames_ = contact.contact_frames_;

        vector_t force_des(contact1.GetNumContacts()*3);
        int j = 0;
        for (int i = 0; i < contact1.in_contact_.size(); i++) {
            if (contact1.in_contact_.at(i)) {
                force_des.segment<3>(3*j) = traj_.GetForce(i, time);
                j++;
            }
        }
        // TODO: Consider removing the logging and/or redoing it
        for (int ee = 0; ee < 4; ee++) {
            log_file_ << std::setw(15) << traj_.GetForce(ee, time).transpose();
        }
        mpc_res_mut_.unlock(); // -------- End Accessing MPC Result -------- //

        vector_t control_action = vector_t::Zero(3*num_inputs_);    // joints, joint velocities, torques

        qp_controller_.UpdateTargetConfig(q_des_);
        qp_controller_.UpdateTargetVel(v_des_);
        qp_controller_.UpdateForceTargets(force_des);

        qp_controller_.UpdateDesiredContacts(contact1);

//        traj_ = traj;

        run_num++;

        log_file_ << std::endl;

        timer.StopTimer();
//        timer.PrintElapsedTime();

        return qp_controller_.ComputeControlAction(q, v, a, contact1, time);
    }

    vector_t MPCController::ReconstructState(const controller::vector_t& q, const controller::vector_t& v,
                                             const controller::vector_t& a) const {

        vector_t state = vector_t::Zero(mpc_.GetModel()->GetNumManifoldStates());

        // Floating base position
        state.head<POS_VARS>() = q.head<3>();

        // COM momentum
        state.segment<POS_VARS>(mpc::SingleRigidBodyModel::LIN_MOM_START) = v.head<3>() * mpc_.GetModel()->GetMass(); // * mpc_.GetModel()->GetMass();

        // Orientation
        Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(q.segment<4>(3)));
        // Note the warning on the pinocchio function!
        pinocchio::quaternion::firstOrderNormalize(quat);
        state(mpc::SingleRigidBodyModel::ORIENTATION_START) = quat.x();
        state(mpc::SingleRigidBodyModel::ORIENTATION_START + 1) = quat.y();
        state(mpc::SingleRigidBodyModel::ORIENTATION_START + 2) = quat.z();
        state(mpc::SingleRigidBodyModel::ORIENTATION_START + 3) = quat.w();

        // Angular Momentum
        Eigen::Quaterniond vel_quat;
        pinocchio::quaternion::exp3(v.segment<POS_VARS>(3), vel_quat);
        Eigen::Quaterniond orientation(static_cast<Eigen::Vector4d>(q.segment<4>(3)));
        vel_quat = orientation.inverse()*vel_quat;
        Eigen::Vector3d vel_frame = pinocchio::quaternion::log3(vel_quat);

        state.segment<POS_VARS>(mpc::SingleRigidBodyModel::ANG_VEL_START + 1) = //.setZero();
                model_.GetIr() * v.segment<POS_VARS>(3); //v.segment<POS_VARS>(3); //model_.GetIr() *

//        std::cout << "vel frame: " << vel_frame.transpose() << std::endl;
//        std::cout << "v: " << v.segment<3>(3).transpose() << std::endl;
//        std::cout << "vcom: " << pin_data_->vcom[0].transpose() << std::endl;
//        std::cout << "I*v: " << (model_.GetIr()*v.segment<3>(3)).transpose() << std::endl;
//        std::cout << "I*vel frame: " << (model_.GetIr()*vel_frame).transpose() << std::endl;
//        std::cout << "I*vcom: " << (model_.GetIr()*pin_data_->vcom[0]).transpose() << std::endl;

//        std::cout << "lin v: " << v.head<3>().transpose() << std::endl;
//        std::cout << "v*m: " << model_.GetMass()*v.head<3>().transpose() << std::endl;

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
        double time = 0;
        std::vector<mpc::vector_3t> ee_locations(4);
        Contact contact;

        int run_num = 0;
        double prev_cost = 1e10;
        double cost_red = 0;

        while(true) {
            state_time_mut_.lock(); // Get time
            double time_act = time_;
            state_time_mut_.unlock();

            if (time_act - time >= info_.integrator_dt) {
                state_time_mut_.lock(); // Get time and state
                state = state_;
                time = time_;
                ee_locations = ee_locations_;
                contact = contact_;
                state_time_mut_.unlock();

                // Run MPC
                mpc_.AdjustForCurrentContacts(time, contact);

                if (!(run_num % gait_opt_freq_) && run_num > 0) {
                    // Print contact times
                    std::cout << "Time: " << time << std::endl;
                    const auto& contact_times = mpc_.GetTrajectory().GetContactTimes();
                    for (int ee = 0; ee < contact_times.size(); ee++) {
                        std::cout << "End effector: " << ee << ".   ";
                        for (int node = 0; node < contact_times.at(ee).size(); node++) {
                            std::cout << std::setw(5) << contact_times.at(ee).at(node).GetTime() << " | ";
                        }
                        std::cout << std::endl;
                    }

                    prev_cost = GaitOptLS(cost_red, time, ee_locations);

                    // Print contact times
                    std::cout << "Time: " << time << std::endl;
                    const auto& contact_times2 = mpc_.GetTrajectory().GetContactTimes();
                    for (int ee = 0; ee < contact_times2.size(); ee++) {
                        std::cout << "End effector: " << ee << ".   ";
                        for (int node = 0; node < contact_times2.at(ee).size(); node++) {
                            std::cout << std::setw(5) << contact_times2.at(ee).at(node).GetTime() << " | ";
                        }
                        std::cout << std::endl;
                    }

                } else {
                    // Print contact times
                    std::cout << "Time: " << time << std::endl;
                    const auto& contact_times = mpc_.GetTrajectory().GetContactTimes();
                    for (int ee = 0; ee < contact_times.size(); ee++) {
                        std::cout << "End effector: " << ee << ".   ";
                        for (int node = 0; node < contact_times.at(ee).size(); node++) {
                            std::cout << std::setw(5) << contact_times.at(ee).at(node).GetTime() << " | ";
                        }
                        std::cout << std::endl;
                    }
                    mpc_.GetRealTimeUpdate(state, time, ee_locations, false);
                }

                // TODO: Look into the contacts - I feel like something is up here

                mpc::Trajectory traj = mpc_.GetTrajectory();
                cost_red = prev_cost - mpc_.GetCost();


                // IK on the trajectory
                utils::Timer ik_timer("traj ik");
                ik_timer.StartTimer();
                FullBodyTrajUpdate(traj);
                ik_timer.StopTimer();
                ik_timer.PrintElapsedTime();

                mpc_res_mut_.lock(); // Set results
                traj_ = traj;
                UpdateTrajViz();
                mpc_res_mut_.unlock();

                run_num++;
            }
        }
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
        const int num_to_update = 10; // traj.GetStates().size();
        for (int i = 0; i < num_to_update; i++) {
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
                vector_t vel(model_.GetFullModelConfigSpace() - 1);
                vel.head<6>() << traj.GetState(i).segment<3>(3)/model_.GetMass(), model_.GetIrInv()*traj.GetState(i).segment<3>(10); //model_.GetIr().inverse()*
                vel.tail(num_inputs_) = (traj.GetFullConfig(i) - traj.GetFullConfig(i-1)).tail(num_inputs_)/info_.integrator_dt;

                traj.UpdateFullVelocity(i, vel);
            }
        }
        traj.UpdateFullVelocity(0, vector_t::Zero(pin_model_.nv));
    }

    void MPCController::GetTargetsFromTraj(const mpc::Trajectory& traj, double time) {
        if (time < traj.GetTime(0)) {
            time = traj.GetTime(0);
        }

        const int nodes_ahead = 1;

        int node = traj.GetNode(time)+nodes_ahead; // TODO: +1 or +0

        q_des_ = traj.GetFullConfig(node);
        v_des_ = traj.GetFullVelocity(node);

        // TODO: Do I want to use the ZoH or the full spline force?
        for (int ee = 0; ee < 4; ee++) {
            force_des_.at(ee) = traj.GetForce(ee, time + nodes_ahead*info_.integrator_dt);
        }
    }

    std::vector<Eigen::Vector2d> MPCController::GetEEBoxCenter() {
        // TODO: Deal with multithreading
        return mpc_.GetEEBoxCenter();
    }

    double MPCController::GaitOptLS(double cost_red, double time, const std::vector<Eigen::Vector3d>& ee_locations) {
        utils::Timer gait_opt_timer("gait opt + derivatives");
        gait_opt_timer.StartTimer();

        const mpc::Trajectory prev_traj = mpc_.GetTrajectory();

        // Create a new MPC and solve with it
        if (mpc_.ComputeDerivativeTerms()) {

            gait_opt_.SetContactTimes(mpc_.GetTrajectory().GetContactTimes());
            gait_opt_.UpdateSizes(mpc_.GetNumDecisionVars(), mpc_.GetNumConstraints());
            double original_cost = mpc_.GetCost();

            utils::Timer qp_partials_timer("qp partials");
            qp_partials_timer.StartTimer();
            mpc_.GetQPPartials(gait_opt_.GetQPPartials());
            qp_partials_timer.StopTimer();
            qp_partials_timer.PrintElapsedTime();


            utils::Timer partials_timer("param partials");
            partials_timer.StartTimer();
            for (int ee = 0; ee < 4; ee++) {
                gait_opt_.SetNumContactTimes(ee, prev_traj.GetNumContactNodes(ee));
                for (int idx = 0; idx < prev_traj.GetNumContactNodes(ee); idx++) {
                    mpc_.ComputeParamPartialsClarabel(prev_traj, gait_opt_.GetParameterPartials(ee, idx), ee, idx);
                }
            }
            partials_timer.StopTimer();
            partials_timer.PrintElapsedTime();

            gait_opt_.ModifyQPPartials(mpc_.GetQPSolution());

            utils::Timer cost_fcn_timer("cost fcn gradient");
            cost_fcn_timer.StartTimer();
            gait_opt_.ComputeCostFcnDerivWrtContactTimes();
            cost_fcn_timer.StopTimer();
            cost_fcn_timer.PrintElapsedTime();

            gait_opt_.OptimizeContactTimes(time, cost_red);

            // Apply the minimizing contact time
            std::vector<mpc::time_v> contact_times;
            double cost_min;
            std::tie(contact_times, cost_min) = gait_opt_.LineSearch(mpc_);
//        mpc.UpdateContactTimes(contact_times);

            gait_opt_timer.StopTimer();
            gait_opt_timer.PrintElapsedTime();

            return cost_min;
        } else {
            std::cerr << "Can't perform gait optimization because MPC was not solved to tolerance." << std::endl;
            gait_opt_timer.StopTimer();
            gait_opt_timer.PrintElapsedTime();
            return 1e30;
        }
    }

} // controller