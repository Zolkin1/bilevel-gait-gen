//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"
#include "mpc.h"
#include "spline/spline.h"
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

double RunGaitOpt(mpc::MPCSingleRigidBody& mpc, mpc::GaitOptimizer& gait_opt, const mpc::Trajectory& prev_traj,
                  double cost_red, double time) {
    double prev_cost = INFINITY;
    if (mpc.ComputeDerivativeTerms()) {

        gait_opt.SetContactTimes(mpc.GetTrajectory().GetContactTimes());
        gait_opt.UpdateSizes(mpc.GetNumDecisionVars(), mpc.GetNumConstraints());
        prev_cost = mpc.GetCost();

        mpc.GetQPPartials(gait_opt.GetQPPartials());
        for (int ee = 0; ee < 4; ee++) {
            gait_opt.SetNumContactTimes(ee, prev_traj.GetNumContactNodes(ee));
            for (int idx = 0; idx < prev_traj.GetNumContactNodes(ee); idx++) {
                mpc.ComputeParamPartialsClarabel(prev_traj, gait_opt.GetParameterPartials(ee, idx), ee, idx);
            }
        }

        gait_opt.ModifyQPPartials(mpc.GetQPSolution());
        gait_opt.ComputeCostFcnDerivWrtContactTimes();

        gait_opt.OptimizeContactTimes(time, cost_red);

        mpc.UpdateContactTimes(gait_opt.GetContactTimes());
    } else {
        std::cerr << "Can't perform gait optimization because MPC was not solved to tolerance." << std::endl;
    }

    return prev_cost;
}

void PrintContactSched(const std::vector<mpc::time_v>& contact_times) {
    for (int ee = 0; ee < contact_times.size(); ee++) {
        std::cout << "End effector: " << ee << ".   ";
        for (int node = 0; node < contact_times.at(ee).size(); node++) {
            std::cout << std::setw(5) << contact_times.at(ee).at(node).GetTime() << " | ";
        }
        std::cout << std::endl;
    }
}

void MPCWithFixedPosition(mpc::MPCSingleRigidBody& mpc, mpc::GaitOptimizer& gait_opt, const vector_t& init_state,
                          std::vector<Eigen::Vector3d> ee_locations, const std::string& robot_xml,
                          const vector_t& standing, const std::unique_ptr<simulator::SimulationRobot>& robot,
                          const mpc::MPCInfo& info, double viz_rate, bool run_gait_opt, simulation::Visualizer viz) {
    mpc.CreateInitialRun(init_state, ee_locations);
    mpc.PrintStats();
    mpc::Trajectory prev_traj = mpc.GetTrajectory();

    vector_t state = standing;
    double prev_cost = mpc.GetCost();
    double cost_red = 0; // TODO: Pick a better number

    const auto contact_sched = mpc.GetTrajectory().GetContactTimes();
    double total_cost = 0;

    const int N = 200;
    for (int i = 0; i < N; i++) {
        double time = i*info.integrator_dt;

        // Get full state through IK
        state = mpc.GetFullTargetState(time, state);
//        vector_t temp = mpc.GetModel()->ConvertTangentStateToManifoldState(prev_traj.GetState(1),
//                                                                           prev_traj.GetState(0));
//        state << temp.head<3>(), temp.segment<4>(6), standing.tail(12);

        // Get end effector locations from trajectory
        for (int j = 0; j < ee_locations.size(); j++) {
            ee_locations.at(j) = prev_traj.GetEndEffectorLocation(j, time);
        }

        // Update Viz
        viz.UpdateState(robot->ConvertPinocchioConfigToMujoco(state)); //mpc.GetTargetConfig(i*info.integrator_dt)));
        viz.GetTrajViz(mpc.CreateVizData(), info.ee_box_size, mpc.GetEEBoxCenter());
        viz.UpdateViz(viz_rate);

        // Gait optimization
        if (run_gait_opt) {
            if (!(i % 10)) {
                prev_cost = RunGaitOpt(mpc, gait_opt, prev_traj, cost_red, time);
//                PrintContactSched(mpc.GetTrajectory().GetContactTimes());
            }
        }

        // Run next MPC
        prev_traj = mpc.GetRealTimeUpdate(prev_traj.GetState(1), time, ee_locations, false);
        cost_red = prev_cost - mpc.GetCost();
        total_cost += mpc.GetCost();
    }

    // Print the final trajectory to a file for viewing
    mpc::Trajectory traj = mpc.GetTrajectory();

    mpc.PrintStats();

    std::cout << "Using gait opt: " << run_gait_opt << std::endl;

    std::cout << "Original contact schedule: " << std::endl;
    PrintContactSched(contact_sched);

    std::cout << "Optimized contact schedule: " << std::endl;
    PrintContactSched(mpc.GetTrajectory().GetContactTimes());

    std::cout << "Average Cost: " << total_cost/static_cast<double>(N) << std::endl;
}

mpc::MPCSingleRigidBody CreateMPC(const mpc::MPCInfo& info, utils::ConfigParser& config, const std::vector<vector_t>& warm_start,
                                  const vector_t& mpc_des_state, const vector_t& init_state,
                                  const std::vector<Eigen::Vector3d>& ee_locations) {
    mpc::MPCSingleRigidBody mpc(info, config.ParseString("robot_urdf"));

    // Set warm starts and defaults
    mpc.SetDefaultGaitTrajectory(mpc::Gaits::Trot, config.ParseNumber<int>("num_polys"), ee_locations);
    mpc.SetStateTrajectoryWarmStart(warm_start);

    matrix_t Q(config.ParseEigenVector("Q_srbd_diag").asDiagonal());

    // Desried state in the lie algebra
    const vector_t des_alg = mpc.GetModel()->ConvertManifoldStateToTangentState(mpc_des_state, init_state);

    // Add in costs
    mpc.AddQuadraticTrackingCost(des_alg, Q);
    mpc.AddForceCost(config.ParseNumber<double>("force_cost"));  // Note: NEED to adjust this based on the number of nodes otherwise it is out-weighed
    mpc.SetQuadraticFinalCost(1*Q);
    mpc.SetLinearFinalCost(-1*Q*des_alg);

    return mpc;
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
    info.real_time_iters = config.ParseNumber<int>("run_time_iterations");
    switch (config.ParseNumber<int>("mpc_verbosity")) {
        case 0:
            info.verbose = mpc::Nothing;
            break;
        case 1:
            info.verbose = mpc::Timing;
            break;
        case 2:
            info.verbose = mpc::Optimization;
            break;
        case 3:
            info.verbose = mpc::All;
            break;
        default:
            throw std::runtime_error("Not a valid verbosity level for MPC.");
    }


    // Read in the inital config and parse it for MPC.
    const vector_t standing = config.ParseEigenVector("init_config");
    const vector_t init_state = config.ParseEigenVector("srb_init");
//    init_state.tail(standing.size()) = standing;

    // Create the warm start
    const std::vector<vector_t> warm_start(info.num_nodes+1, init_state);

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

    std::vector<Eigen::Vector3d> ee_locations(4);
    for (int i = 0; i < ee_locations.size(); i++) {
        for (int j = 0; j < 3; j++) {
            ee_locations.at(i)(j) = ee_pos.at(i).at(j);
        }
    }

    mpc::MPCSingleRigidBody mpc1 = CreateMPC(info, config, warm_start, mpc_des_state, init_state, ee_locations);
    mpc::MPCSingleRigidBody mpc2 = CreateMPC(info, config, warm_start, mpc_des_state, init_state, ee_locations);

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
                                                                 config.ParseEigenVector("Q_srbd_diag").asDiagonal());

    // Make the robot for visualization
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);


    mpc::GaitOptimizer gait_optimizer1(4, 10, 10, 10, 1, 0.05);
    mpc::GaitOptimizer gait_optimizer2(4, 10, 10, 10, 1, 0.05);

    // Visualize results
    simulation::Visualizer viz(robot_file);
    robot->SetSimModel(viz.GetModel());

    // MPC w/ fixed position
    MPCWithFixedPosition(mpc1, gait_optimizer1, init_state, ee_locations, config.ParseString("robot_xml"), standing,
                         robot, info, config.ParseNumber<double>("viz_rate"), false, viz);


    MPCWithFixedPosition(mpc2, gait_optimizer2, init_state, ee_locations, config.ParseString("robot_xml"), standing,
                         robot, info, config.ParseNumber<double>("viz_rate"), true, viz);

    // MPC w/ fixed position + line search

    // MPC w/ high precision solves


}