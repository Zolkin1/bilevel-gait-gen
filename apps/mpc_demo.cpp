//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "config_parser.h"
#include "mpc.h"

using vector_t = Eigen::VectorXd;

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

    mpc::MPC mpc(info, config.ParseString("robot_urdf"));

    vector_t state(6+6+12);
    state << 1, 0, 0.1, 0, 0, 0, 0, 0.35, 1, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    mpc.Solve(state);
}