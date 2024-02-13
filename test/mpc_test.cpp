//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>

#include "config_parser.h"
#include "spline/spline.h"
#include "models/centroidal_model.h"
#include "mpc_single_rigid_body.h"

#include "qp/clarabel_interface.h"

TEST_CASE("Basic MPC", "[mpc]") {

    std::string config_file("/home/zolkin/AmberLab/bilevel-gait-gen/test/a1_configuration_test.yaml");
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
    info.nom_state = config.ParseEigenVector("init_config");

    mpc::MPCSingleRigidBody mpc(info, config.ParseString("robot_urdf"));

}

TEST_CASE("Transformations", "[mpc][utils]") {

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

    Eigen::VectorXd state(25);
    Eigen::VectorXd ref_state(25);
    state << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0.1968, 0.1958, 0.0811, 0.9573, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25;
    ref_state << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 1, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;

    for (int i = 0; i < state.size(); i++) {
        REQUIRE_THAT(state(i) - mpc::CentroidalModel::ConvertAlgebraStateToManifoldState(
                mpc::CentroidalModel::ConvertManifoldStateToAlgebraState(state, ref_state), ref_state)(i),
                     WithinAbs(0, MARGIN));
    }
}

TEST_CASE("Normal Spline", "[mpc][spline]") {
    using namespace mpc;
    using Catch::Matchers::WithinAbs;

    std::vector<double> times;
    times.push_back(0.35);
    times.push_back(0.75);

    constexpr double MARGIN = 5e-3;

    // ----------- Basic Creation ----------- //
    const int num_polys1 = 2;
    Spline spline1 = Spline(num_polys1, times, true, Spline::Normal);

    REQUIRE(spline1.GetPolyTimes().size() == 4);
    REQUIRE(spline1.GetEndTime() == 0.75);
    REQUIRE(spline1.GetNumConstant() == 2);
    REQUIRE(spline1.GetTotalPolyVars() == 4);
    REQUIRE(spline1.GetNumPolyTimes() == 4);

    const int num_polys2 = 3;
    Spline spline2 = Spline(num_polys2, times, true, Spline::Normal);

    REQUIRE(spline2.GetPolyTimes().size() == 5);
    REQUIRE(spline2.GetEndTime() == 0.75);
    REQUIRE(spline2.GetNumConstant() == 2);
    REQUIRE(spline2.GetTotalPolyVars() == 6);
    REQUIRE(spline2.GetNumPolyTimes() == 5);

    SECTION("Assigning values") {
        const double position = 1;
        spline1.SetAllPositions(position);
        for (int i = 0; i < spline1.GetPolyVars().size(); i++) {
            REQUIRE(spline1.GetPolyVars().at(i).at(0) == position);
        }

        std::vector<double> vars = {3,4};
        spline1.SetPolyVars(1, vars);
        REQUIRE(spline1.GetPolyVars().at(1).at(0) == vars.at(0));
        REQUIRE(spline1.GetPolyVars().at(1).at(1) == vars.at(1));


        // Set a constant
        std::vector<double> vars2 = {5};
        spline1.SetPolyVars(2, vars2);
        REQUIRE(spline1.GetPolyVars().at(2).at(0) == vars2.at(0));
        REQUIRE(spline1.GetPolyVars().at(3).at(0) == vars2.at(0));
    }

    SECTION("Checking values") {
        std::vector<double> end1 = {1, 4};
        spline1.SetPolyVars(1, end1);

        REQUIRE(spline1.ValueAt(0) == 0);
        REQUIRE_THAT(spline1.ValueAt(0.175), WithinAbs(1, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.030172), WithinAbs(0.06170, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.144827), WithinAbs(0.83841, MARGIN));
        REQUIRE(spline1.ValueAt(0.6) == 0);

        std::vector<double> end2 = {2, -3};
        spline2.SetPolyVars(1, end1);
        spline2.SetPolyVars(2, end2);
        REQUIRE_THAT(spline2.ValueAt(0.11666), WithinAbs(1, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.23333), WithinAbs(2, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.036 + 0.11666666), WithinAbs(1.32076, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.076 + 0.11666666), WithinAbs(1.85302, MARGIN + 0.05));   // Slightly more numerical error
        REQUIRE_THAT(spline2.ValueAt(0.112 + 0.11666666), WithinAbs(2.00823, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.35), WithinAbs(0.0, MARGIN));

        std::vector<double> end3 = {5};
        spline2.SetPolyVars(3, end3);
        REQUIRE_THAT(spline2.ValueAt(0.35), WithinAbs(5, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.55), WithinAbs(5, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.75), WithinAbs(5, MARGIN));
    }

    SECTION("Linearization/Coefficients") {
        std::vector<double> end1 = {1, 4};
        std::vector<double> end2 = {2, -3};
        std::vector<double> end3 = {5};
        spline2.SetPolyVars(1, end1);
        spline2.SetPolyVars(2, end2);
        spline2.SetPolyVars(3, end3);

        vector_t vars(6);
        vars << 0, 1, 4, 2, -3, 5;

        double time = 0;
        while (time <= spline2.GetEndTime()) {
            vector_t vars_coef = spline2.GetPolyVarsLin(time);

            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(time);
            REQUIRE(vars_affecting == vars_coef.size());

            double lin_val = vars_coef.dot(vars.segment(vars_index-vars_affecting,vars_affecting));
            REQUIRE_THAT(spline2.ValueAt(time), WithinAbs(lin_val, MARGIN));
            time += 0.015;
        }


        std::vector<double> end4 = {3, -0.1};
        std::vector<double> end5 = {2.1};
        spline1.SetPolyVars(0, end5);
        spline1.SetPolyVars(1, end4);

        vector_t vars2(4);
        vars2 << 2.1, 3, -0.1, 0;

        time = 0;
        while (time <= spline1.GetEndTime()) {
            vector_t vars_coef = spline1.GetPolyVarsLin(time);

            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = spline1.GetVarsIndexEnd(time);
            REQUIRE(vars_affecting == vars_coef.size());

            double lin_val = vars_coef.dot(vars2.segment(vars_index-vars_affecting,vars_affecting));
            REQUIRE_THAT(spline1.ValueAt(time), WithinAbs(lin_val, MARGIN));
            time += 0.015;
        }

    }

    SECTION("Adding and removing spline parts") {
        std::vector<double> end1 = {1, 4};
        std::vector<double> end2 = {2, -3};
        std::vector<double> end3 = {5};
        spline2.SetPolyVars(1, end1);
        spline2.SetPolyVars(2, end2);
        spline2.SetPolyVars(3, end3);

        REQUIRE(spline2.GetPolyTimes().size() == 5);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 6);
        REQUIRE(spline2.GetNumPolyTimes() == 5);

        spline2.RemoveUnused(0.1);
        REQUIRE(spline2.GetPolyTimes().size() == 5);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 6);
        REQUIRE(spline2.GetNumPolyTimes() == 5);

        spline2.RemoveUnused(0.2);
        REQUIRE(spline2.GetPolyTimes().size() == 4);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 1);
        REQUIRE(spline2.GetTotalPolyVars() == 5);
        REQUIRE(spline2.GetNumPolyTimes() == 4);

        int vars_index, vars_affecting;
        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(0.2);
        REQUIRE(vars_index == 4);
        REQUIRE(vars_affecting == 4);

        spline2.RemoveUnused(0.4);
        REQUIRE(spline2.GetPolyTimes().size() == 2);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 1);
        REQUIRE(spline2.GetTotalPolyVars() == 1);
        REQUIRE(spline2.GetNumPolyTimes() == 2);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(0.71);
        REQUIRE(vars_index == 1);
        REQUIRE(vars_affecting == 1);

        spline2.RemoveUnused(0.41);
        REQUIRE(spline2.GetPolyTimes().size() == 2);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 1);
        REQUIRE(spline2.GetTotalPolyVars() == 1);
        REQUIRE(spline2.GetNumPolyTimes() == 2);

        spline2.AddPoly(0.5);
        REQUIRE(spline2.GetPolyTimes().size() == 5);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(1.25, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 6);
        REQUIRE(spline2.GetNumPolyTimes() == 5);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(0.71);
        REQUIRE(vars_index == 1);
        REQUIRE(vars_affecting == 1);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(0.9);
        REQUIRE(vars_index == 3);
        REQUIRE(vars_affecting == 3);

        spline2.AddPoly(0.5);
        REQUIRE(spline2.GetPolyTimes().size() == 6);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(1.75, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 6);
        REQUIRE(spline2.GetNumPolyTimes() == 6);

        spline2.AddPoly(0.5);
        REQUIRE(spline2.GetPolyTimes().size() == 9);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(2.25, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 3);
        REQUIRE(spline2.GetTotalPolyVars() == 11);
        REQUIRE(spline2.GetNumPolyTimes() == 9);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(1);
        REQUIRE(vars_index == 5);
        REQUIRE(vars_affecting == 4);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(1.1);
        REQUIRE(vars_index == 6);
        REQUIRE(vars_affecting == 3);

        spline2.AddPoly(0.35);
        REQUIRE(spline2.GetPolyTimes().size() == 10);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(2.6, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 3);
        REQUIRE(spline2.GetTotalPolyVars() == 11);
        REQUIRE(spline2.GetNumPolyTimes() == 10);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(1.5);
        REQUIRE(vars_index == 6);
        REQUIRE(vars_affecting == 1);

        REQUIRE(spline2.ValueAt(2.6) == 0);
        REQUIRE(spline2.ValueAt(0.6) == 5);
    }
}

TEST_CASE("Force Spline", "[mpc][spline]") {
    using namespace mpc;
    using Catch::Matchers::WithinAbs;

    std::vector<double> times;
    times.push_back(0.35);
    times.push_back(0.75);

    constexpr double MARGIN = 5e-3;

    // ----------- Basic Creation ----------- //
    const int num_polys1 = 2;
    Spline spline1 = Spline(num_polys1, times, true, Spline::Force);

    REQUIRE(spline1.GetPolyTimes().size() == 4);
    REQUIRE(spline1.GetEndTime() == 0.75);
    REQUIRE(spline1.GetNumConstant() == 2);
    REQUIRE(spline1.GetTotalPolyVars() == 2);
    REQUIRE(spline1.GetNumPolyTimes() == 4);

    const int num_polys2 = 3;
    Spline spline2 = Spline(num_polys2, times, true, Spline::Force);

    REQUIRE(spline2.GetPolyTimes().size() == 5);
    REQUIRE(spline2.GetEndTime() == 0.75);
    REQUIRE(spline2.GetNumConstant() == 2);
    REQUIRE(spline2.GetTotalPolyVars() == 4);
    REQUIRE(spline2.GetNumPolyTimes() == 5);

    SECTION("Assigning values") {
        const double position = 1;
        spline1.SetAllPositions(position);
        for (int i = 0; i < spline1.GetPolyVars().size(); i++) {
            if (spline1.GetPolyVars().at(i).size() == 2) {
                REQUIRE(spline1.GetPolyVars().at(i).at(0) == position);
            } else {
                REQUIRE(spline1.GetPolyVars().at(i).at(0) == 0);
            }
        }

        std::vector<double> vars = {3,4};
        spline1.SetPolyVars(1, vars);
        REQUIRE(spline1.GetPolyVars().at(1).at(0) == vars.at(0));
        REQUIRE(spline1.GetPolyVars().at(1).at(1) == vars.at(1));
    }

    SECTION("Checking values") {
        std::vector<double> end1 = {1, 4};
        spline1.SetPolyVars(1, end1);

        REQUIRE(spline1.ValueAt(0) == 0);
        REQUIRE_THAT(spline1.ValueAt(0.175), WithinAbs(1, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.030172), WithinAbs(0.06170, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.144827), WithinAbs(0.83841, MARGIN));
        REQUIRE(spline1.ValueAt(0.6) == 0);

        std::vector<double> end2 = {2, -3};
        spline2.SetPolyVars(1, end1);
        spline2.SetPolyVars(2, end2);
        REQUIRE_THAT(spline2.ValueAt(0.11666), WithinAbs(1, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.23333), WithinAbs(2, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.036 + 0.11666666), WithinAbs(1.32076, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.076 + 0.11666666), WithinAbs(1.85302, MARGIN + 0.05));   // Slightly more numerical error
        REQUIRE_THAT(spline2.ValueAt(0.112 + 0.11666666), WithinAbs(2.00823, MARGIN));
        REQUIRE_THAT(spline2.ValueAt(0.35), WithinAbs(0.0, MARGIN));
    }

    SECTION("Linearization/Coefficients") {
        std::vector<double> end1 = {1, 4};
        std::vector<double> end2 = {2, -3};
        spline2.SetPolyVars(1, end1);
        spline2.SetPolyVars(2, end2);

        vector_t vars(4);
        vars << 1, 4, 2, -3;

        double time = 0;
        while (time <= spline2.GetEndTime()) {
            if (!spline2.IsConstant(time)) {
                vector_t vars_coef = spline2.GetPolyVarsLin(time);

                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(time);
                REQUIRE(vars_affecting == vars_coef.size());

                double lin_val = vars_coef.dot(vars.segment(vars_index-vars_affecting,vars_affecting));
                REQUIRE_THAT(spline2.ValueAt(time), WithinAbs(lin_val, MARGIN));
            }
            time += 0.015;
        }


        std::vector<double> end4 = {3, -0.1};
        spline1.SetPolyVars(1, end4);

        vector_t vars2(2);
        vars2 << 3, -0.1;

        time = 0;
        while (time <= spline1.GetEndTime()) {
            if (!spline1.IsConstant(time)) {
                vector_t vars_coef = spline1.GetPolyVarsLin(time);

                int vars_index, vars_affecting;
                std::tie(vars_index, vars_affecting) = spline1.GetVarsIndexEnd(time);
                REQUIRE(vars_affecting == vars_coef.size());

                double lin_val = vars_coef.dot(vars2.segment(vars_index - vars_affecting, vars_affecting));
                REQUIRE_THAT(spline1.ValueAt(time), WithinAbs(lin_val, MARGIN));
            }
            time += 0.015;
        }
    }

    SECTION("Adding and removing spline parts") {
        std::vector<double> end1 = {1, 4};
        std::vector<double> end2 = {2, -3};
        spline2.SetPolyVars(1, end1);
        spline2.SetPolyVars(2, end2);

        REQUIRE(spline2.GetPolyTimes().size() == 5);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 4);
        REQUIRE(spline2.GetNumPolyTimes() == 5);

        spline2.RemoveUnused(0.1);
        REQUIRE(spline2.GetPolyTimes().size() == 5);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 4);
        REQUIRE(spline2.GetNumPolyTimes() == 5);

        spline2.RemoveUnused(0.2);
        REQUIRE(spline2.GetPolyTimes().size() == 4);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 1);
        REQUIRE(spline2.GetTotalPolyVars() == 4);
        REQUIRE(spline2.GetNumPolyTimes() == 4);

        int vars_index, vars_affecting;
        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(0.2);
        REQUIRE(vars_index == 4);
        REQUIRE(vars_affecting == 4);

        spline2.RemoveUnused(0.4);
        REQUIRE(spline2.GetPolyTimes().size() == 2);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 1);
        REQUIRE(spline2.GetTotalPolyVars() == 0);
        REQUIRE(spline2.GetNumPolyTimes() == 2);

        spline2.RemoveUnused(0.41);
        REQUIRE(spline2.GetPolyTimes().size() == 2);
        REQUIRE(spline2.GetEndTime() == 0.75);
        REQUIRE(spline2.GetNumConstant() == 1);
        REQUIRE(spline2.GetTotalPolyVars() == 0);
        REQUIRE(spline2.GetNumPolyTimes() == 2);

        spline2.AddPoly(0.5);
        REQUIRE(spline2.GetPolyTimes().size() == 5);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(1.25, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 4);
        REQUIRE(spline2.GetNumPolyTimes() == 5);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(0.9);
        REQUIRE(vars_index == 2);
        REQUIRE(vars_affecting == 2);

        spline2.AddPoly(0.5);
        REQUIRE(spline2.GetPolyTimes().size() == 6);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(1.75, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 2);
        REQUIRE(spline2.GetTotalPolyVars() == 4);
        REQUIRE(spline2.GetNumPolyTimes() == 6);

        spline2.AddPoly(0.5);
        REQUIRE(spline2.GetPolyTimes().size() == 9);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(2.25, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 3);
        REQUIRE(spline2.GetTotalPolyVars() == 8);
        REQUIRE(spline2.GetNumPolyTimes() == 9);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(1);
        REQUIRE(vars_index == 4);
        REQUIRE(vars_affecting == 4);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(1.1);
        REQUIRE(vars_index == 4);
        REQUIRE(vars_affecting == 2);

        spline2.AddPoly(0.35);
        REQUIRE(spline2.GetPolyTimes().size() == 10);
        REQUIRE_THAT(spline2.GetEndTime(), WithinAbs(2.6, MARGIN));
        REQUIRE(spline2.GetNumConstant() == 3);
        REQUIRE(spline2.GetTotalPolyVars() == 8);
        REQUIRE(spline2.GetNumPolyTimes() == 10);

        std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(1.15);
        REQUIRE(vars_index == 4);
        REQUIRE(vars_affecting == 2);

        REQUIRE(spline2.ValueAt(2.6) == 0);
        spline2.SetPolyVars(2, {2,2});
        REQUIRE_THAT(spline2.ValueAt(0.916), WithinAbs(2.0, MARGIN));
    }
}

TEST_CASE("Constant Splines", "[mpc][spline]") {
    using namespace mpc;
    using Catch::Matchers::WithinAbs;

    std::vector<double> times;
    times.push_back(0.35);
    times.push_back(0.75);

    constexpr double MARGIN = 5e-3;

    // ----------- Basic Creation ----------- //
    const int num_polys1 = 2;
    Spline spline1 = Spline(num_polys1, times, true, Spline::Constants);

    Spline spline2 = Spline(num_polys1, times, false, Spline::Constants);

    REQUIRE(spline1.GetPolyTimes().size() == 3);
    REQUIRE(spline1.GetEndTime() == 0.75);
    REQUIRE(spline1.GetNumConstant() == 3);
    REQUIRE(spline1.GetTotalPolyVars() == 2);
    REQUIRE(spline1.GetNumPolyTimes() == 3);

    for (int i = 0; i < spline1.GetPolyTimes().size(); i++) {
        REQUIRE(spline1.GetPolyVars().at(i).size() == 1);
    }


    SECTION("Assigning values") {
        const double position = 1;
        spline1.SetAllPositions(position);
        for (int i = 0; i < spline1.GetPolyVars().size(); i++) {
            REQUIRE(spline1.GetPolyVars().at(i).at(0) == 1);
        }

        std::vector<double> vars = {3};
        spline1.SetPolyVars(1, vars);
        REQUIRE(spline1.GetPolyVars().at(1).at(0) == vars.at(0));
    }

    SECTION("Checking values") {
        std::vector<double> end1 = {4};
        spline1.UpdatePolyVar(1, end1);

        REQUIRE(spline1.ValueAt(0) == 0);
        REQUIRE_THAT(spline1.ValueAt(0.0965517), WithinAbs(0.745254, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.217241), WithinAbs(2.71007, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.35), WithinAbs(4, MARGIN));
        REQUIRE_THAT(spline1.ValueAt(0.75), WithinAbs(0.0, MARGIN));
    }

    SECTION("Linearization/Coefficients") {
        std::vector<double> end1 = {1};
        std::vector<double> end2 = {2};
        spline1.UpdatePolyVar(0, end1);
        spline1.UpdatePolyVar(1, end2);
        spline1.UpdatePolyVar(2, end2);

        vector_t vars(2);
        vars << 1, 2;

        double time = 0;
        while (time <= spline1.GetEndTime()) {
            vector_t vars_coef = spline1.GetPolyVarsLin(time);

            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = spline1.GetVarsIndexEnd(time);
            REQUIRE(vars_affecting == vars_coef.size());

            double lin_val = vars_coef.dot(vars.segment(vars_index-vars_affecting,vars_affecting));
            REQUIRE_THAT(spline1.ValueAt(time), WithinAbs(lin_val, MARGIN));

            time += 0.015;
        }

        spline2.UpdatePolyVar(0, end1);
        spline2.UpdatePolyVar(1, end1);
        spline2.UpdatePolyVar(2, end2);
        time = 0;
        while (time <= spline2.GetEndTime()) {
            vector_t vars_coef = spline2.GetPolyVarsLin(time);

            int vars_index, vars_affecting;
            std::tie(vars_index, vars_affecting) = spline2.GetVarsIndexEnd(time);
            REQUIRE(vars_affecting == vars_coef.size());

            double lin_val = vars_coef.dot(vars.segment(vars_index-vars_affecting,vars_affecting));
            REQUIRE_THAT(spline2.ValueAt(time), WithinAbs(lin_val, MARGIN));

            time += 0.015;
        }
    }
}

void GetData(mpc::QPData& data) {
    data.num_decision_vars = 3;
    data.num_dynamics_constraints = 2;
    if (data.using_clarabel_) {
        data.num_force_box_constraints_ = 4;
    } else {
        data.num_force_box_constraints_ = 2;
    }

    data.InitQPMats();

    mpc::vector_t P_vals(3);
    P_vals << 3.001, 4, 0.5;

    mpc::vector_t q(3);
    q << 0.1, 4.6, 2;

    data.cost_mat_.SetVectorDiagonally(P_vals, 0, 0);
    data.cost_linear = q;

    mpc::matrix_t A(2,3);
    A << 1, 1, 0, 1.3, 0, 0.2;

    mpc::vector_t b(2);
    b << 1, 3;

    data.constraint_mat_.SetMatrix(A, 0, 0);
    data.dynamics_constants = b;

    A << -2, 0, 0.9, 1, 8, 5;
    b << 3.1, 13.3;

    mpc::vector_t b_lb(2);
    b_lb << -2, -5;

    data.constraint_mat_.SetMatrix(A, 2, 0);
    if (data.using_clarabel_) {
        data.constraint_mat_.SetMatrix(-A, 4, 0);
    }
    data.force_box_lb_ = b_lb;
    data.force_box_ub_ = b;

    data.ConstructVectors();
    data.ConstructSparseMats();

    std::cout << "Constraint mat: \n" << data.sparse_constraint_.toDense() << std::endl;
    std::cout << "upper bound: \n" << data.ub_ << std::endl;
}

void CheckVectorsEqual(const mpc::vector_t& vec1, const mpc::vector_t& vec2, const double MARGIN) {
    using Catch::Matchers::WithinAbs;

    REQUIRE(vec1.size() == vec2.size());

    for (int i = 0; i < vec1.size(); i++) {
        REQUIRE_THAT(vec1(i) - vec2(i), WithinAbs(0, MARGIN));
    }
}

TEST_CASE("Clarabel Solver", "[mpc][clarabel]") {

    constexpr double MARGIN = 1e-4;

    using namespace mpc;
    using Catch::Matchers::WithinAbs;

    QPData clar_data(false, 10, 10, {Dynamics, ForceBox});
    clar_data.using_clarabel_ = true;

    QPData osqp_data(false, 10, 10, {Dynamics, ForceBox});
    osqp_data.using_clarabel_ = false;

    GetData(clar_data);
    GetData(osqp_data);

    ClarabelInterface clarabel(clar_data, true);
    clarabel.SetupQP(clar_data, {});
    vector_t clara_sol = clarabel.Solve(clar_data);

    vector_t clara_dx = clarabel.Computedx(clar_data.sparse_cost_, clar_data.cost_linear, clara_sol);
    vector_t zero = vector_t::Zero(1);
    clarabel.SetupDerivativeCalcs(clara_dx, zero, zero, clar_data);

    sp_matrix_t dP(clar_data.sparse_cost_.rows(), clar_data.sparse_cost_.cols());
    sp_matrix_t dA(clar_data.num_dynamics_constraints, 3);
    sp_matrix_t dG(clar_data.num_force_box_constraints_, 3);
    clarabel.CalcDerivativeWrtMats(dP, dA, dG);


    OSQPInterface osqp(osqp_data, true);
    vector_t zero_warmstart = vector_t::Zero(3);
    osqp.SetupQP(osqp_data, zero_warmstart);
    vector_t osqp_sol = osqp.Solve(osqp_data);

    REQUIRE(osqp.GetSolveQuality() == clarabel.GetSolveQuality());

    CheckVectorsEqual(osqp_sol, clara_sol, MARGIN);

    osqp.Computedx(osqp_data.sparse_cost_, osqp_data.cost_linear, osqp_sol);
    vector_t osqp_dx = osqp.Getdx();

    CheckVectorsEqual(osqp_dx, clara_dx, MARGIN);

    vector_t zero_b = vector_t::Zero(osqp_data.GetTotalNumConstraints());
    osqp.SetupDerivativeCalcs(osqp_dx, zero_b, zero_b);

    sp_matrix_t osqp_dP(clar_data.sparse_cost_.rows(), clar_data.sparse_cost_.cols());
    sp_matrix_t osqp_dA(2, 3);
    osqp.CalcDerivativeWrtMats(osqp_dP, osqp_dA);

    std::cout << "clara dual: " << clarabel.GetDualSolution().transpose() << std::endl;
    std::cout << "osqp dual: " << osqp.GetDualSolution().transpose() << std::endl;

    matrix_t c_dP = dP.toDense();
    matrix_t o_dP = osqp_dP.toDense();
    std::cout << "Clara dP: " << c_dP << std::endl;
    std::cout << "OSQP dP: " << o_dP << std::endl;

    for (int i = 0; i < clar_data.num_decision_vars; i++) {
        for (int j = 0; j < clar_data.num_decision_vars; j++) {
//            REQUIRE_THAT(c_dP(i, j) - o_dP(i, j), WithinAbs(0, MARGIN));
        }
    }

    matrix_t c_dA = dA.toDense();
    matrix_t c_dG = dG.toDense();
    matrix_t o_dA = osqp_dA.toDense();

    std::cout << "Clara dA: \n" << c_dA << std::endl;
    std::cout << "Clara dG: \n" << c_dG << std::endl;
    std::cout << "OSQP dA: \n" << o_dA << std::endl;

    // TODO: OSQP dA only matches mine where the sparsity pattern is non-zero! (my implementation is correct afaik)

    for (int i = 0; i < osqp_data.num_dynamics_constraints; i++) {
        for (int j = 0; j < osqp_data.num_decision_vars; j++) {
            if (osqp_data.sparse_constraint_.toDense()(i,j) != 0) {   // Due to OSQP bug can only expect this to match when the term is non-zero
                REQUIRE_THAT(c_dA(i, j) - o_dA(i, j), WithinAbs(0, MARGIN));
            }
        }
    }

    for (int i = 0; i < osqp_data.num_force_box_constraints_; i++) {
        for (int j = 0; j < osqp_data.num_decision_vars; j++) {
            int idx = i + osqp_data.num_dynamics_constraints;
            if (osqp_data.sparse_constraint_.toDense()(idx, j) != 0) { // Due to OSQP bug can only expect this to match when the term is non-zero
                REQUIRE_THAT(c_dG(i, j) + c_dG(i + osqp_data.num_force_box_constraints_, j)
                             - o_dA(idx, j), WithinAbs(0, MARGIN));
            }
        }
    }

}