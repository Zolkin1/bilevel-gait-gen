//
// Created by zolkin on 1/30/24.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>
#include <iostream>

#include "end_effector_splines.h"

TEST_CASE("End Effector Splines", "[splines]") {
    constexpr double MARGIN = 1e-3;

    using namespace mpc;
    using Catch::Matchers::WithinAbs;

    std::vector<double> times;
    const int num_contacts = 10;
    const bool start_constant = true;
    const int num_force_polys = 3;

    for (int i = 0; i < num_contacts; i++) {
        times.emplace_back(0.2*i);
    }

    std::vector<mpc::EndEffectorSplines> splines;
    splines.emplace_back(num_contacts, times, !start_constant, num_force_polys);
    splines.emplace_back(num_contacts, times, start_constant, num_force_polys);


    SECTION("Setting Vars") {
        for (auto& spline : splines) {
            for (int coord = 0; coord < 3; coord++) {
                const std::vector<double> spline_times = spline.GetTimes();

                const std::vector<int> position_nodes = spline.GetMutableNodes(EndEffectorSplines::Position);
                for (auto& it: position_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 2;
                    spline.SetVars(EndEffectorSplines::Position, coord, it, vars);
                    REQUIRE(spline.ValueAt(EndEffectorSplines::Position, coord, spline_times.at(it)) == it);
                }

                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 2;
                    spline.SetVars(EndEffectorSplines::Force, coord, it, vars);
                    REQUIRE(spline.ValueAt(EndEffectorSplines::Force, coord, spline_times.at(it)) == it);
                }
            }
        }
    }

    SECTION("Checking Values") {

        for (int coord = 0; coord < 3; coord++) {
            const std::vector<int> position_nodes = splines.at(0).GetMutableNodes(EndEffectorSplines::Position);
            for (auto& it: position_nodes) {
                Eigen::Vector2d vars;
                vars << it, it - 1;
                splines.at(0).SetVars(EndEffectorSplines::Position, coord, it, vars);
            }
            REQUIRE(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0) == 0);
            REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0.103448), WithinAbs(0.525852, MARGIN));
            REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0.503448), WithinAbs(3.10341, MARGIN));

            const std::vector<int> position_nodes2 = splines.at(1).GetMutableNodes(EndEffectorSplines::Position);
            for (auto& it: position_nodes2) {
                Eigen::Vector2d vars;
                vars << it, it - 1;
                splines.at(1).SetVars(EndEffectorSplines::Position, coord, it, vars);
            }
            REQUIRE(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0) == 0);
            REQUIRE_THAT(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0.103448), WithinAbs(0.0, MARGIN));
            REQUIRE_THAT(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0.25517), WithinAbs(0.74525, MARGIN));

            const std::vector<int> force_nodes = splines.at(0).GetMutableNodes(EndEffectorSplines::Force);
            for (auto& it: force_nodes) {
                Eigen::Vector2d vars;
                vars << it, it - 1;
                splines.at(0).SetVars(EndEffectorSplines::Force, coord, it, vars);
            }
            REQUIRE(splines.at(0).ValueAt(EndEffectorSplines::Force, coord, 0) == 0);
            REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Force, coord, 0.103448), WithinAbs(0.0, MARGIN));
            REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Force, coord, 0.26666 + 0.0551724),
                         WithinAbs(2.90697, MARGIN));
        }
    }

    SECTION("Linearization/Coefficients") {
        for (auto& spline : splines) {
            for (int coord = 0; coord < 3; coord++) {
                const std::vector<double> spline_times = spline.GetTimes();
                const double total_time = spline.GetEndTime();

                const std::vector<int> position_nodes = spline.GetMutableNodes(EndEffectorSplines::Position);
                for (auto& it: position_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 3.1;
                    spline.SetVars(EndEffectorSplines::Position, coord, it, vars);
                }

                const double N = 100.0;
                const vector_t spline_as_vec = spline.GetSplineAsQPVec(EndEffectorSplines::Position, coord);
                for (int i = 0; i < N; i++) {
                    const double time = static_cast<double>(i) * (total_time / N);
                    int vars_idx, vars_effecting;
                    std::tie(vars_idx, vars_effecting) = spline.GetVarsIdx(EndEffectorSplines::Position, coord, time);
                    const vector_t vars_lin = spline.GetPolyVarsLin(EndEffectorSplines::Position, coord, time);
                    REQUIRE(vars_lin.size() == vars_effecting);
                    const double val = spline_as_vec.segment(vars_idx, vars_effecting).dot(vars_lin);
                    REQUIRE_THAT(spline.ValueAt(EndEffectorSplines::Position, coord, time), WithinAbs(val, MARGIN));
                }

                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 1.4;
                    spline.SetVars(EndEffectorSplines::Force, coord, it, vars);
                }

                const vector_t spline_as_vec2 = spline.GetSplineAsQPVec(EndEffectorSplines::Force, coord);
                for (int i = 0; i < N; i++) {
                    const double time = static_cast<double>(i) * (total_time / N);

                    if (spline.IsForceMutable(time)) {
                        int vars_idx, vars_effecting;
                        std::tie(vars_idx, vars_effecting) = spline.GetVarsIdx(EndEffectorSplines::Force, coord, time);
                        const vector_t vars_lin = spline.GetPolyVarsLin(EndEffectorSplines::Force, coord, time);
                        REQUIRE(vars_lin.size() == vars_effecting);
                        const double val = spline_as_vec2.segment(vars_idx, vars_effecting).dot(vars_lin);
                        REQUIRE_THAT(spline.ValueAt(EndEffectorSplines::Force, coord, time), WithinAbs(val, MARGIN));
                    } else {
                        REQUIRE(spline.ValueAt(EndEffectorSplines::Force, coord, time) == 0);
                    }
                }
            }
        }
    }

    SECTION("Adding and removing polys") {
        const double addt_time = 0.2;
        for (auto& spline : splines) {
            // ---------------- Adding ---------------- //
            for (int i = 0; i < 3; i++) {
                spline.AddPoly(addt_time);
                REQUIRE_THAT(spline.GetEndTime(),
                             WithinAbs(times.at(times.size() - 1) + (i + 1) * addt_time, MARGIN));
            }

            for (int coord = 0; coord < 3; coord++) {
                // Checking values after adding polys
                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << it - 1, .75;
                    spline.SetVars(EndEffectorSplines::Force, coord, it, vars);
                }

                const double total_time = spline.GetEndTime();
                const double N = 100.0;
                const vector_t spline_as_vec2 = spline.GetSplineAsQPVec(EndEffectorSplines::Force, coord);
                for (int i = 0; i < N; i++) {
                    const double time = static_cast<double>(i) * (total_time / N);

                    if (spline.IsForceMutable(time)) {
                        int vars_idx, vars_effecting;
                        std::tie(vars_idx, vars_effecting) = spline.GetVarsIdx(EndEffectorSplines::Force, coord, time);
                        const vector_t vars_lin = spline.GetPolyVarsLin(EndEffectorSplines::Force, coord, time);
                        REQUIRE(vars_lin.size() == vars_effecting);
                        const double val = spline_as_vec2.segment(vars_idx, vars_effecting).dot(vars_lin);
                        REQUIRE_THAT(spline.ValueAt(EndEffectorSplines::Force, coord, time), WithinAbs(val, MARGIN));
                    } else {
                        REQUIRE(spline.ValueAt(EndEffectorSplines::Force, coord, time) == 0);
                    }
                }
            }
        }

        // ---------------- Removing ---------------- //
        for (auto& spline : splines) {
            spline.RemovePoly(0.5);
            REQUIRE_THAT(spline.GetEndTime(), WithinAbs(times.at(times.size() - 1) + (3) * addt_time, MARGIN));

            for (int coord = 0; coord < 3; coord++) {

                // Now check values again
                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << 2 * it - 1, .5;
                    spline.SetVars(EndEffectorSplines::Force, coord, it, vars);
                }

                const double total_time = spline.GetEndTime();
                const double N = 100.0;
                const vector_t spline_as_vec2 = spline.GetSplineAsQPVec(EndEffectorSplines::Force, coord);
                for (int i = 0; i < N; i++) {
                    const double time =
                            static_cast<double>(i) * ((total_time - spline.GetStartTime()) / N) + spline.GetStartTime();

                    if (spline.IsForceMutable(time)) {
                        int vars_idx, vars_effecting;
                        std::tie(vars_idx, vars_effecting) = spline.GetVarsIdx(EndEffectorSplines::Force, coord, time);
                        const vector_t vars_lin = spline.GetPolyVarsLin(EndEffectorSplines::Force, coord, time);
                        REQUIRE(vars_lin.size() == vars_effecting);
                        const double val = spline_as_vec2.segment(vars_idx, vars_effecting).dot(vars_lin);
                        REQUIRE_THAT(spline.ValueAt(EndEffectorSplines::Force, coord, time), WithinAbs(val, MARGIN));
                    } else {
                        REQUIRE(spline.ValueAt(EndEffectorSplines::Force, coord, time) == 0);
                    }
                }
            }
        }
    }
}