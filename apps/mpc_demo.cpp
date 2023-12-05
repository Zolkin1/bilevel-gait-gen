//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"
#include "mpc.h"
#include "inputs.h"
#include "spline.h"
#include "mpc_controller.h"

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
    info.time_horizon = config.ParseNumber<double>("time_horizon");
    info.num_qp_iterations = config.ParseNumber<int>("num_qp");
    info.friction_coef = config.ParseNumber<double>("friction_coef");
    info.vel_bounds = config.ParseEigenVector("vel_bounds");
    info.joint_bounds = config.ParseEigenVector("joint_bounds");
    info.ee_frames = config.ParseStdVector<std::string>("collision_frames");
    info.num_switches = config.ParseNumber<int>("num_switches");
    info.integrator_dt = config.ParseNumber<double>("integrator_dt");
    info.num_contacts = info.ee_frames.size();


    vector_t standing = config.ParseEigenVector("standing_config");
    standing.segment<4>(3) = ConvertMujocoQuatToPinocchioQuat(standing.segment<4>(3));
    vector_t state = vector_t::Zero(6 + standing.size());
    state.tail(standing.size()) = standing;

    std::vector<vector_t> warm_start_states;
    for (int i = 0; i < info.num_nodes+1; i++) {
        warm_start_states.push_back(state);
    }

    mpc::MPC mpc(info, config.ParseString("robot_urdf"));

    auto switching_times = mpc::MPC::CreateDefaultSwitchingTimes(info.num_switches, 4, info.time_horizon);
    mpc::Trajectory traj(info.num_nodes+1, config.ParseEigenVector("init_config").size() + 6, 12,
                         switching_times,
                         info.time_horizon/info.num_nodes);

    vector_t standing_state = state;
    for (int i = 0; i < info.num_nodes+1; i++) {
        state(1) += i/2.0;
        traj.SetState(i, state);
    }

    mpc::Inputs input(switching_times, config.ParseEigenVector("init_config").size() + 6 - 13,
                      info.num_nodes, info.time_horizon/info.num_nodes);

    // Make the same force for each ee to ease of testing right now
    std::array<mpc::Spline, 3> forces1 = {mpc::Spline(3, switching_times.at(0), true),
                                         mpc::Spline(3, switching_times.at(0), true),
                                         mpc::Spline(3, switching_times.at(0), true)};
    std::array<mpc::Spline, 3> forces2 = {mpc::Spline(3, switching_times.at(0), true),
                                          mpc::Spline(3, switching_times.at(0), true),
                                          mpc::Spline(3, switching_times.at(0), true)};
    std::array<mpc::Spline, 3> positions1 = {mpc::Spline(3, switching_times.at(0), false),
                                            mpc::Spline(3, switching_times.at(0), false),
                                         mpc::Spline(3, switching_times.at(0), false)};
    std::array<mpc::Spline, 3> positions2 = {mpc::Spline(3, switching_times.at(0), false),
                                            mpc::Spline(3, switching_times.at(0), false),
                                            mpc::Spline(3, switching_times.at(0), false)};
    // For the test just make everything constant, non-zero

    for (int ee = 0; ee < 4; ee++) {
        traj.SetEndEffectorSplines(ee, forces1.at(0), positions1.at(0));

//        for (int coord = 0; coord < 3; coord++) {
//            for (int poly = 0; poly < input.GetForces().at(ee).at(coord).GetNumPolyTimes(); poly++) {
//                std::vector<double> vars(input.GetForces().at(ee).at(coord).GetNumPolyVars(poly));
//                vars.at(0) = 1;
//                if (vars.size() == 2) {
//                    vars.at(1) = 0; // Constants for now
//                }
//                if (ee%2 == 0) {
//                    forces1.at(coord).SetPolyVars(poly, vars);
//                } else {
//                    forces2.at(coord).SetPolyVars(poly, vars);
//                }
//            }
//        for (int poly = 0; poly < input.GetPositions().at(0).at(0).GetTotalPoly(); poly++) {
//            positions.at(coord).SetPolyVars(poly, vars);
//        }
//        }
    }

    for (int ee = 0; ee < 4; ee++) {
        if (ee%2 == 0) {
            input.SetEndEffectorForce(ee, forces1);
        } else {
            input.SetEndEffectorForce(ee, forces2);
        }
//        input.SetEndEffectorPosition(ee, positions);
    }

    traj.SetInput(input);

    std::vector<std::array<double, 3>> ee_pos;
    ee_pos.push_back({0.2, 0.2, 0});
    ee_pos.push_back({0.2, -0.2, 0});
    ee_pos.push_back({-0.2, 0.2, 0});
    ee_pos.push_back({-0.2, 0.2, 0});

    for (int ee = 0; ee < 4; ee++) {
        traj.SetPositionsForAllTime(ee, ee_pos.at(ee));
    }

    mpc.SetWarmStartTrajectory(traj);

    // TODO: work on speed
    vector_t curr_state(6+7+12);
    curr_state = standing_state;
    curr_state(7) += 0.5;
    //    curr_state << 1, 0, 0.1,
//            0, 0, 0,
//            0, 0, 0.35,
//            1, 0, 0, 0,
//            0, 0, 1, 0, 0.1, 0, 0, 0, 0.2, 0, 0, 0;

    vector_t state_des = vector_t::Zero(6+7+12);
//    state_des(11) = 1;
    state_des = standing_state;
    state_des(7) = 1;
    state_des(6) = 1;

    matrix_t Q = matrix_t::Identity(24, 24);
    Q.topLeftCorner<6,6>() = matrix_t::Zero(6,6);
    Q(6,6) = 30;
    Q(7,7) = 30;
    Q(8,8) = 10;

    const vector_t des_alg = mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(state_des, standing_state);
    std::cout << des_alg << std::endl;

//    vector_t w = vector_t::Zero(24);
    mpc.AddQuadraticTrackingCost(des_alg, Q);
//    mpc.SetQuadraticFinalCost(10*Q);
//    mpc.SetLinearFinalCost(-10*Q*des_alg);

//    curr_state.segment<4>(9) = ConvertMujocoQuatToPinocchioQuat(curr_state.segment<4>(9));
    mpc::Trajectory solve_traj = mpc.Solve(curr_state, 0);
    solve_traj.PrintTrajectoryToFile("solved_traj.txt");


    for (int i = 0; i < 3; i++) {
        mpc.Solve(curr_state, 0);//(i+1)*info.integrator_dt);
    }
}