//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>

#include "mpc.h"
#include "config_parser.h"

TEST_CASE("basic mpc", "[mpc]") {

    std::string config_file("a1_configuration_test.yaml");
    utils::ConfigParser config = utils::ConfigParser(config_file);

    mpc::MPCInfo info;
    info.discretization_steps = config.ParseNumber<int>("discretization_steps");
    info.num_nodes = config.ParseNumber<int>("num_nodes");
    info.time_horizon = config.ParseNumber<int>("time_horizon");
    info.num_qp_iterations = config.ParseNumber<int>("num_qp");
    info.friction_coef = config.ParseNumber<double>("friction_coef");
    info.vel_bounds = config.ParseEigenVector("vel_bounds");
    info.joint_bounds = config.ParseEigenVector("joint_bounds");
    info.ee_frames = config.ParseStdVector<std::string>("collision_frames");

    mpc::MPC mpc(info, config.ParseString("robot_urdf"));

}

TEST_CASE("transformations", "[mpc][utils]") {

    using Catch::Matchers::WithinAbs;

    double constexpr MARGIN = 1e-3;

    Eigen::Vector4d quat;
    quat << 0.7071, 0, 0, 0.7071;
    Eigen::Vector3d rot;
    rot << 0, 0, 1.57078;

    Eigen::Vector3d ZYXRotSol = mpc::CentroidalModel::ConvertQuaternionToZYXRot(quat);
    for (int i = 0; i < ZYXRotSol.size(); i++) {
        REQUIRE_THAT(ZYXRotSol(i), WithinAbs(rot(i), MARGIN));
    }

    quat << 0.36515, 0.54772, 0.7303, 0.18257;
    rot << 2.3562, -0.3398, 1.4289;
    ZYXRotSol = mpc::CentroidalModel::ConvertQuaternionToZYXRot(quat);
    for (int i = 0; i < ZYXRotSol.size(); i++) {
        REQUIRE_THAT(ZYXRotSol(i), WithinAbs(rot(i), MARGIN));
    }

    quat << 0.5773, 0.5773, 0, 0.5773;
    rot << 1.1069, 0.72957, 2.03423;
    ZYXRotSol = mpc::CentroidalModel::ConvertQuaternionToZYXRot(quat);
    for (int i = 0; i < ZYXRotSol.size(); i++) {
        REQUIRE_THAT(ZYXRotSol(i), WithinAbs(rot(i), MARGIN));
    }

    Eigen::Vector4d quat_sol = mpc::CentroidalModel::ConvertZYXRotToQuaternion(rot);
    for (int i = 0; i < quat_sol.size(); i++) {
        REQUIRE_THAT(quat_sol(i), WithinAbs(quat(i), MARGIN));
    }

    rot << 0.25, 0.35, 0.45;
    quat << 0.1968, 0.1958, 0.0811, 0.9573;
    quat_sol = mpc::CentroidalModel::ConvertZYXRotToQuaternion(rot);
    for (int i = 0; i < quat_sol.size(); i++) {
        REQUIRE_THAT(quat_sol(i), WithinAbs(quat(i), MARGIN));
    }
}

// TODO: Write a spline test
