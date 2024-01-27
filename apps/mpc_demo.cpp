//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"
#include "mpc.h"
#include "spline.h"
#include "mpc_controller.h"
#include "visualization.h"
#include "simulation_robot.h"
#include "gait_optimizer.h"

using vector_t = Eigen::VectorXd;
using matrix_t = Eigen::MatrixXd;

Eigen::Vector4d ConvertMujocoQuatToPinocchioQuat(const Eigen::Vector4d& quat) {
    Eigen::Vector4d pin_quat;
    pin_quat(3) = quat(0);
    pin_quat(0) = quat(1);
    pin_quat(1) = quat(2);
    pin_quat(2) = quat(3);

    return pin_quat;
}

int main() {
    // TODO: when I make the file local via the cmake file then I only get it copied over after a compilation
    std::string config_file("/home/zolkin/AmberLab/bilevel-gait-gen/apps/a1_configuration.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    mpc::MPCInfo info;
    info.discretization_steps = config.ParseNumber<double>("discretization_steps");
    info.num_nodes = config.ParseNumber<int>("num_nodes");
    info.num_qp_iterations = config.ParseNumber<int>("num_qp");
    info.friction_coef = config.ParseNumber<double>("friction_coef");
    info.vel_bounds = config.ParseEigenVector("vel_bounds");
    info.joint_bounds_lb = config.ParseEigenVector("joint_bounds_lb");
    info.joint_bounds_ub = config.ParseEigenVector("joint_bounds_ub");
    info.ee_frames = config.ParseStdVector<std::string>("collision_frames");
    info.num_switches = config.ParseNumber<int>("num_switches");
    info.integrator_dt = config.ParseNumber<double>("integrator_dt");
    info.num_contacts = info.ee_frames.size();
    info.force_bound = config.ParseNumber<double>("force_bound");
    info.swing_height = config.ParseNumber<double>("swing_height");
    info.foot_offset = config.ParseNumber<double>("foot_offset");
    info.nom_state = config.ParseEigenVector("init_config");
    info.ee_box_size = config.ParseEigenVector("ee_box_size");

    mpc::MPCSingleRigidBody mpc(info, config.ParseString("robot_urdf"));


    // Read in the inital config and parse it for MPC.
    vector_t standing = config.ParseEigenVector("init_config");
    vector_t init_state = config.ParseEigenVector("srbd_init");
//    init_state.tail(standing.size()) = standing;

    // Create the warm start
    std::vector<vector_t> warm_start(info.num_nodes+1, init_state);

    // Create the goal state
    vector_t mpc_des_state = init_state;
    mpc_des_state.head<2>() << config.ParseNumber<double>("x_des"), config.ParseNumber<double>("y_des");
    mpc_des_state.segment<2>(3) << config.ParseNumber<double>("xdot_des"), config.ParseNumber<double>("ydot_des");

    // Inital guess end effector positions
    std::array<std::array<double, 3>, 4> ee_pos{};
    ee_pos.at(0) = {0.2, 0.2, 0};
    ee_pos.at(1) = {0.2, -0.2, 0};
    ee_pos.at(2) = {-0.2, 0.2, 0};
    ee_pos.at(3) = {-0.2, -0.2, 0};

    // Set warm starts and defaults
    mpc.SetDefaultGaitTrajectory(mpc::Gaits::Trot, config.ParseNumber<int>("num_polys"), ee_pos);
    mpc.SetStateTrajectoryWarmStart(warm_start);

    matrix_t Q(config.ParseEigenVector("Q_srbd_diag").asDiagonal());

    std::cout << "Q: " << Q << std::endl;

    // Desried state in the lie algebra
    const vector_t des_alg = mpc.GetModel()->ConvertManifoldStateToTangentState(mpc_des_state, init_state);
            //mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(mpc_des_state, init_state);
    std::cout << des_alg << std::endl;

    // Add in costs
    mpc.AddQuadraticTrackingCost(des_alg, Q);
    mpc.AddForceCost(config.ParseNumber<double>("force_cost"));  // Note: NEED to adjust this based on the number of nodes otherwise it is out-weighed
    mpc.SetQuadraticFinalCost(1*Q);
    mpc.SetLinearFinalCost(-1*Q*des_alg);

    // Create the MPC controller (only used here for the visualizer)
    std::unique_ptr<controller::Controller> mpc_controller;
    mpc_controller = std::make_unique<controller::MPCController>(config.ParseNumber<double>("control_rate"),
                                                                 config.ParseString("robot_urdf"),
                                                                 config.ParseString("foot_type"),
                                                                 config.ParseEigenVector("init_vel").size(),
                                                                 config.ParseEigenVector("torque_bounds"),
                                                                 config.ParseNumber<double>("friction_coef"),
                                                                 config.ParseStdVector<double>("base_pos_gains"),
                                                                 config.ParseStdVector<double>("base_ang_gains"),
                                                                 config.ParseEigenVector("kp_joint_gains"),
                                                                 config.ParseEigenVector("kd_joint_gains"),
                                                                 config.ParseNumber<double>("leg_tracking_weight"),
                                                                 config.ParseNumber<double>("torso_tracking_weight"),
                                                                 config.ParseNumber<double>("force_tracking_weight"),
                                                                 info,
                                                                 warm_start,
                                                                 mpc_des_state,
                                                                 config.ParseNumber<int>("num_polys"),
                                                                 Q);

    // Make the robot for visualization
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);

    std::vector<Eigen::Vector3d> ee_locations(4);
    for (int i = 0; i < ee_locations.size(); i++) {
        for (int j = 0; j < 3; j++) {
            ee_locations.at(i)(j) = ee_pos.at(i).at(j);
        }
    }

    mpc.CreateInitialRun(init_state, ee_locations);
    mpc.PrintStats();
    mpc::Trajectory prev_traj = mpc.GetTrajectory();

    mpc::GaitOptimizer gait_optimizer(4, 10, 10, 10, 1, 0.05);

    // Visualize results
    simulation::Visualizer viz(config.ParseString("robot_xml"));
    robot->SetSimModel(viz.GetModel());
    vector_t state = standing;
    for (int i = 0; i < info.num_nodes + 200; i++) {
        // Only grabbing the COM states for now (for debugging)
//        state.head<7>() = mpc.GetFullTargetState(i*info.integrator_dt, state).head<7>();
        state = mpc.GetFullTargetState(i*info.integrator_dt, state);
        for (int j = 0; j < ee_locations.size(); j++) {
            ee_locations.at(j) = prev_traj.GetEndEffectorLocation(j, (i)*info.integrator_dt);
        }

        viz.UpdateState(robot->ConvertPinocchioConfigToMujoco(state)); //mpc.GetTargetConfig(i*info.integrator_dt)));
        viz.GetTrajViz(mpc.CreateVizData(), info.ee_box_size, mpc.GetEEBoxCenter());
        viz.UpdateViz(config.ParseNumber<double>("viz_rate"));
//        mpc.GetRealTimeUpdate(config.ParseNumber<int>("run_time_iterations"),
//                prev_traj.GetState(1), i*info.integrator_dt, ee_locations);
        prev_traj = mpc.GetRealTimeUpdate(config.ParseNumber<int>("run_time_iterations"),
                                          prev_traj.GetState(1), i*info.integrator_dt, ee_locations);
        // Gait optimization
//        if (i == 0) {
//            gait_optimizer.SetContactTimes(mpc.GetTrajectory().GetContactTimes());
//        }
        // TODO: Is this being performed with linearizations from the new trajectory and/or is that messing up the calcs?
//        gait_optimizer.UpdateSizes(mpc.GetNumDecisionVars(), mpc.GetNumConstraints());
//        if (mpc.ComputeDerivativeTerms()) {
//            mpc.GetQPPartials(gait_optimizer.GetPartials());
//            const mpc::Trajectory traj = mpc.GetTrajectory();
//            for (int ee = 0; ee < 4; ee++) {
//                gait_optimizer.SetNumContactTimes(ee, traj.GetNumContactNodes(ee));
//                for (int idx = 0; idx < traj.GetNumContactNodes(ee); idx++) {
//                    mpc.ComputeParamPartials(prev_traj, gait_optimizer.GetParameterPartials(ee, idx), ee, idx);
//                }
//            }
//
//            gait_optimizer.ModifyQPPartials(mpc.GetQPSolution());
//            gait_optimizer.ComputeCostFcnDerivWrtContactTimes();
//
//            gait_optimizer.OptimizeContactTimes();

//                mpc.UpdateContactTimes(gait_optimizer_.GetContactTimes());

//        mpc.GetRealTimeUpdate(config.ParseNumber<int>("run_time_iterations"), temp_state, i*info.integrator_dt);
//        }

//        prev_traj = mpc.GetTrajectory();
    }

    // Print the final trajectory to a file for viewing
    mpc::Trajectory traj = mpc.GetTrajectory();
    traj.PrintTrajectoryToFile("demo_final_traj.txt");

    mpc.PrintStats();

}