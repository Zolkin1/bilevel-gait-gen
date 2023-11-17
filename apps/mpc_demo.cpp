//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"
#include "mpc.h"
#include "inputs.h"
#include "spline.h"

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
    std::string config_file("a1_configuration.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    mpc::MPCInfo info;
    info.discretization_steps = config.ParseNumber<double>("discretization_steps");
    info.num_nodes = config.ParseNumber<int>("num_nodes");
    info.time_horizon = config.ParseNumber<int>("time_horizon");
    info.num_qp_iterations = config.ParseNumber<int>("num_qp");
    info.friction_coef = config.ParseNumber<double>("friction_coef");
    info.vel_bounds = config.ParseEigenVector("vel_bounds");
    info.joint_bounds = config.ParseEigenVector("joint_bounds");
    info.ee_frames = config.ParseStdVector<std::string>("collision_frames");
    info.num_switches = config.ParseNumber<int>("num_switches");
    info.integrator_dt = config.ParseNumber<double>("integrator_dt");
    info.num_contacts = info.ee_frames.size();

    mpc::MPC mpc(info, config.ParseString("robot_urdf"));

    auto switching_times = mpc::MPC::CreateDefaultSwitchingTimes(info.num_switches, 4, info.time_horizon);
    mpc::Trajectory traj(info.num_nodes, config.ParseEigenVector("init_config").size() + 6, 12,
                         switching_times,
                         info.time_horizon/info.num_nodes);

    vector_t standing = config.ParseEigenVector("standing_config");
    standing.segment<4>(3) = ConvertMujocoQuatToPinocchioQuat(standing.segment<4>(3));
    vector_t state = vector_t::Zero(6 + standing.size());
    state.tail(standing.size()) = standing;
    for (int i = 0; i < info.num_nodes; i++) {
        state(1) += i/2.0;
        traj.SetState(i, state);
    }

    mpc::Inputs input(switching_times, config.ParseEigenVector("init_config").size() + 6 - 13,
                      info.num_nodes, info.time_horizon/info.num_nodes);

    // Make the same force for each ee to ease of testing right now
    std::array<mpc::Spline, 3> forces = {mpc::Spline(2, switching_times.at(0), true), mpc::Spline(2, switching_times.at(0), true),
                                         mpc::Spline(2, switching_times.at(0), true)};
    std::array<mpc::Spline, 3> positions = {mpc::Spline(2, switching_times.at(0), false), mpc::Spline(2, switching_times.at(0), false),
                                         mpc::Spline(2, switching_times.at(0), false)};

    // For the test just make everything constant, non-zero
    for (int ee = 0; ee < 4; ee++) {
        for (int coord = 0; coord < 3; coord++) {
            for (int poly = 0; poly < input.GetForces().at(ee).at(coord).GetNumPolyTimes(); poly++) {
                std::vector<double> vars(input.GetForces().at(ee).at(coord).GetNumPolyVars(poly));
                vars.at(0) = 1;
                if (vars.size() == 2) {
                    vars.at(1) = 0; // Constants for now
                }
                forces.at(coord).SetPolyVars(poly, vars);
            }
//        for (int poly = 0; poly < input.GetPositions().at(0).at(0).GetTotalPoly(); poly++) {
//            positions.at(coord).SetPolyVars(poly, vars);
//        }
        }
    }

    // TODO: might want to consider setting feet positions to something that is not all zeros (still zero in z direction)

    for (int ee = 0; ee < 4; ee++) {
        input.SetEndEffectorForce(ee, forces);
//        input.SetEndEffectorPosition(ee, positions);
    }

    traj.SetInput(input);

    mpc.SetWarmStartTrajectory(traj);

    matrix_t Q = matrix_t::Identity(24, 24);
    vector_t w = vector_t::Zero(24);

    mpc.SetLinearCostTerm(w);
    mpc.SetQuadraticCostTerm(Q);
    mpc.SetQuadraticFinalCost(10*Q);

    vector_t curr_state(6+7+12);
    curr_state << 1, 0, 0.1,
                    0, 0, 0,
                    0, 0, 0.35,
                    1, 0, 0, 0,
                    0, 0, 1, 0, 0.1, 0, 0, 0, 0.2, 0, 0, 0;
    curr_state.segment<4>(9) = ConvertMujocoQuatToPinocchioQuat(curr_state.segment<4>(9));
    mpc.Solve(curr_state);
}