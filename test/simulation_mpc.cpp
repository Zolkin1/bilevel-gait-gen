//
// Created by zolkin on 2/19/24.
//

#include "config_parser.h"
#include "mpc.h"
#include "spline/spline.h"
#include "mpc_controller.h"
#include "simple_simulation.h"
#include "simulation_robot.h"
#include "gait_optimizer.h"

using namespace mpc;

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

void PrintContactSched(const std::vector<mpc::time_v>& contact_times) {
    for (int ee = 0; ee < contact_times.size(); ee++) {
        std::cout << "End effector: " << ee << ".   ";
        for (int node = 0; node < contact_times.at(ee).size(); node++) {
            std::cout << std::setw(5) << contact_times.at(ee).at(node).GetTime() << " | ";
        }
        std::cout << std::endl;
    }
}

int main() {
    // TODO: when I make the file local via the cmake file then I only get it copied over after a compilation
    std::string config_file("/home/zolkin/AmberLab/bilevel-gait-gen/apps/a1_configuration.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    const double viz_rate = config.ParseNumber<double>("viz_rate");

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

//    mpc::MPCSingleRigidBody mpc = CreateMPC(info, config, warm_start, mpc_des_state, init_state, ee_locations);

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
                                                                 config.ParseEigenVector("Q_srbd_diag").asDiagonal(),
                                                                 config.ParseNumber<int>("gait_opt_freq"));

    // Make the robot for visualization
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);


//    mpc::GaitOptimizer gait_opt(4, 10, 10, 10, 1, 0.05);

    if (viz_rate != 1.0/config.ParseNumber<double>("control_rate")) {
        throw std::runtime_error("different viz rate and integration dt is not allowed");
    }

    // Set the robot's initial condition
    vector_t init_vel = config.ParseEigenVector("init_vel");
    robot->SetInitialCondition(standing, init_vel);


    // Set up controller solver
    robot->InitController(standing, init_state);
    robot->DefineContacts(config.ParseStringVector("collision_frames"),
                          config.ParseStdVector<int>("collision_bodies"));

    // Visualize results
    simulation::SimpleSimulation sim(robot_file);
    robot->SetSimModel(sim.GetModelPointer());

//    mpc.CreateInitialRun(init_state, ee_locations);
//    mpc::Trajectory prev_traj = mpc.GetTrajectory();

    mjData* data = sim.GetDataPointer();
    for (int i = 0; i < standing.size(); i++) {
        data->qpos[i] = robot->ConvertPinocchioConfigToMujoco(standing)(i);
    }
    for (int i = 0; i < init_vel.size(); i++) {
        data->qvel[i] = robot->ConvertPinocchioVelToMujoco(init_vel)(i);
    }

    // Update Sim
    sim.GetTrajViz(robot->GetTrajViz(), info.ee_box_size, robot->GetEEBoxCenter());
    sim.UpdateSim(viz_rate);

    sim.PauseSim();

    int constexpr N = 6500;
    utils::Timer timer("sim loop");
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++) {
        double time = i*info.integrator_dt;

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        const std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        if (time_span.count() >= viz_rate) {
            t1 = std::chrono::high_resolution_clock::now();

            timer.StartTimer();

            robot->GetControlAction(sim.GetDataPointer(), sim.GetControlPointer());

            // Update Sim
            sim.GetTrajViz(robot->GetTrajViz(), info.ee_box_size, robot->GetEEBoxCenter());
            sim.UpdateSim(viz_rate);
            timer.StopTimer();
        } else {
            const double time_left = -time_span.count() + viz_rate;
            const int time_left_us = static_cast<int>(time_left*1e6);
            std::this_thread::sleep_for(std::chrono::microseconds(time_left_us));
        }

//        std::this_thread::sleep_for(std::chrono::milliseconds (10));

//        timer.PrintElapsedTime();
    }

    // Print the final trajectory to a file for viewing
//    mpc::Trajectory traj = mpc.GetTrajectory();
//
//    mpc.PrintStats();
//
//    std::cout << "Optimized contact schedule: " << std::endl;
//    PrintContactSched(mpc.GetTrajectory().GetContactTimes());
//
//    std::cout << "Average Cost: " << total_cost/static_cast<double>(N) << std::endl;

}