//
// Created by zolkin on 2/19/24.
//

#include "config_parser.h"
#include "mpc.h"
#include "spline/spline.h"
#include "mpc_controller.h"
#include "visualization.h"
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

double GaitOptLS(mpc::MPCSingleRigidBody& mpc, mpc::GaitOptimizer& gait_opt,
                 double cost_red, double time, const MPCInfo& info, utils::ConfigParser& config,
                 const vector_t& mpc_des_state, const std::vector<Eigen::Vector3d>& ee_locations,
                 const std::vector<vector_t>& warm_start, bool fixed_pos) {
    utils::Timer gait_opt_timer("gait opt + derivatives");
    gait_opt_timer.StartTimer();

    double prev_cost = INFINITY;
    const Trajectory prev_traj = mpc.GetTrajectory();

    // Create a new MPC and solve with it
    if (mpc.ComputeDerivativeTerms()) {

        gait_opt.SetContactTimes(mpc.GetTrajectory().GetContactTimes());
        gait_opt.UpdateSizes(mpc.GetNumDecisionVars(), mpc.GetNumConstraints());
        double original_cost = mpc.GetCost();

        utils::Timer qp_partials_timer("qp partials");
        qp_partials_timer.StartTimer();
        mpc.GetQPPartials(gait_opt.GetQPPartials());
        qp_partials_timer.StopTimer();
        qp_partials_timer.PrintElapsedTime();


        utils::Timer partials_timer("param partials");
        partials_timer.StartTimer();
        for (int ee = 0; ee < 4; ee++) {
            gait_opt.SetNumContactTimes(ee, prev_traj.GetNumContactNodes(ee));
            for (int idx = 0; idx < prev_traj.GetNumContactNodes(ee); idx++) {
                mpc.ComputeParamPartialsClarabel(prev_traj, gait_opt.GetParameterPartials(ee, idx), ee, idx);
            }
        }
        partials_timer.StopTimer();
        partials_timer.PrintElapsedTime();

        gait_opt.ModifyQPPartials(mpc.GetQPSolution());

        utils::Timer cost_fcn_timer("cost fcn gradient");
        cost_fcn_timer.StartTimer();
        gait_opt.ComputeCostFcnDerivWrtContactTimes();
        cost_fcn_timer.StopTimer();
        cost_fcn_timer.PrintElapsedTime();

        gait_opt.OptimizeContactTimes(time, cost_red);

        // Apply the minimizing contact time
        std::vector<time_v> contact_times;
        double cost_min;
        std::tie(contact_times, cost_min) = gait_opt.LineSearch(mpc, time, ee_locations, prev_traj.GetState(1));
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

Trajectory MPCTest(mpc::MPCSingleRigidBody& mpc, double time, const std::vector<Eigen::Vector3d>& ee_locations) {
    const Trajectory traj = mpc.GetTrajectory();
    return mpc.GetRealTimeUpdate(traj.GetState(1), time, ee_locations, false);
}

void MPCLineSearch(mpc::MPCSingleRigidBody& mpc, utils::ConfigParser& config, const std::vector<vector_t>& warm_start,
                   const vector_t& mpc_des_state, mpc::GaitOptimizer& gait_opt, const vector_t& init_state,
                   std::vector<Eigen::Vector3d>& ee_locations, const std::string& robot_xml,
                   const vector_t& standing, const std::unique_ptr<simulator::SimulationRobot>& robot,
                   const mpc::MPCInfo& info, double viz_rate, simulation::Visualizer viz,
                   bool fixed_pos) {

    mpc.CreateInitialRun(init_state, ee_locations);
    mpc.PrintStats();
    mpc::Trajectory prev_traj = mpc.GetTrajectory();

    vector_t state = standing;

    double total_cost = 0;
    double prev_cost = 0;
    double cost_red = 0;

    const int N = 300;
    for (int i = 0; i < N; i++) {
        double time = i*info.integrator_dt;
        if (fixed_pos) {
            time = 0;
        }

        // Get end effector locations from trajectory
        for (int j = 0; j < ee_locations.size(); j++) {
            ee_locations.at(j) = prev_traj.GetEndEffectorLocation(j, time);
        }


        // Update Viz
        // Get full state through IK
        state = mpc.GetFullTargetState(time, state);
        viz.UpdateState(robot->ConvertPinocchioConfigToMujoco(state));
        viz.GetTrajViz(mpc.CreateVizData(), info.ee_box_size, mpc.GetEEBoxCenter());
        viz.UpdateViz(viz_rate);

        // TODO: Remove if I don't want the state perturbation
//        if (time == 1) {
//            vector_t state_new = prev_traj.GetState(1);
//            state_new(3) += 6;
//            state_new(4) += 12;
//            prev_traj.SetState(1, state_new);
//            mpc.SetWarmStartTrajectory(prev_traj);
//            std::cout << "State perturbation applied." << std::endl;
//        }

        // Gait optimization
        if (!(i % 50000) && i > 0) { //5
            prev_cost = GaitOptLS(mpc, gait_opt, cost_red, time, info, config, mpc_des_state,
                                  ee_locations, warm_start, fixed_pos);

            prev_traj = mpc.GetTrajectory();
//            prev_traj = mpc.GetRealTimeUpdate(prev_traj.GetState(1), time, ee_locations, false);

        std::cout << "---------" << std::endl;

//            std::cout << "cost diff: " << prev_cost - mpc.GetCost() << std::endl;

 //            std::cout << "Time: " << time << std::endl;
 //            PrintContactSched(mpc.GetTrajectory().GetContactTimes());
        } else {
            prev_traj = mpc.GetRealTimeUpdate(prev_traj.GetState(1), time, ee_locations, false);
        }

        // Run next MPC
//        if (fixed_pos) {
//            prev_traj = mpc.GetRealTimeUpdate(prev_traj.GetState(0), time, ee_locations, false);
//        } else {
//            prev_traj = mpc.GetRealTimeUpdate(prev_traj.GetState(1), time, ee_locations, false);
//        }

        cost_red = prev_cost - mpc.GetCost();
//        std::cout << "cost error (should be 0 if gait optimized): " << cost_red << std::endl;
        total_cost += mpc.GetCost();
    }

    // Print the final trajectory to a file for viewing
    mpc::Trajectory traj = mpc.GetTrajectory();

    mpc.PrintStats();

    std::cout << "Optimized contact schedule: " << std::endl;
    PrintContactSched(mpc.GetTrajectory().GetContactTimes());

    std::cout << "Average Cost: " << total_cost/static_cast<double>(N) << std::endl;
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
    info.force_cost = config.ParseNumber<double>("force_cost");
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

    // Create the goal state
    vector_t mpc_des_state = config.ParseEigenVector("srb_target"); //init_state;
//    mpc_des_state.head<2>() << config.ParseNumber<double>("x_des"), config.ParseNumber<double>("y_des");
//    mpc_des_state.segment<2>(3) << config.ParseNumber<double>("xdot_des"), config.ParseNumber<double>("ydot_des");
//    mpc_des_state.segment<4>(6) << 0, 0, 0, 1;

    // Create the warm start
    std::vector<vector_t> warm_start(info.num_nodes+1); //, init_state);
    for (int i = 0; i < warm_start.size(); i++) {
        warm_start.at(i) = (i/warm_start.size())*(mpc_des_state - init_state) + init_state;
    }

    // Inital guess end effector positions
    std::array<std::array<double, 3>, 4> ee_pos{};
    ee_pos.at(0) = {0.1526, 0.12523, 0.011089};
    ee_pos.at(1) = {0.1526, -0.12523, 0.011089};
    ee_pos.at(2) = {-0.208321844, 0.1363286, 0.01444};
    ee_pos.at(3) = {-0.208321844, -0.1363286, 0.01444};

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
                                                                 config.ParseEigenVector("Q_srbd_diag").asDiagonal(),
                                                                 config.ParseNumber<int>("gait_opt_freq"),
                                                                 config.ParseString("log_file"));

    // Make the robot for visualization
    auto robot_file = config.ParseString("robot_xml");
    std::unique_ptr<simulator::SimulationRobot> robot = std::make_unique<simulator::SimulationRobot>(robot_file, mpc_controller);


    mpc::GaitOptimizer gait_optimizer1(4, 10, 10, 10, 1, 0.05);
    mpc::GaitOptimizer gait_optimizer2(4, 10, 10, 10, 1, 0.05);

    // Visualize results
    simulation::Visualizer viz(robot_file);
    robot->SetSimModel(viz.GetModel());

    // MPC w/ fixed position
    MPCLineSearch(mpc1, config, warm_start, mpc_des_state, gait_optimizer1, init_state, ee_locations, config.ParseString("robot_xml"), standing,
                         robot, info, config.ParseNumber<double>("viz_rate"), viz, false);
}