//
// Created by zolkin on 1/30/24.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>
#include <iostream>

#include "spline/end_effector_splines.h"

TEST_CASE("End Effector Splines", "[splines]") {
    constexpr double MARGIN = 1e-3;

    using namespace mpc;
    using Catch::Matchers::WithinAbs;

    std::vector<double> times;
    const int num_contacts = 5; // 10
    const bool start_constant = true;
    const int num_force_polys = 3;

    for (int i = 0; i < num_contacts; i++) {
        times.emplace_back(0.2*i);
    }

    std::vector<mpc::EndEffectorSplines> splines;
    splines.emplace_back(num_contacts, times, !start_constant, num_force_polys);
    splines.emplace_back(num_contacts, times, start_constant, num_force_polys);

    constexpr double FORCE_MULT = 100;

    SECTION("Setting Vars") {
        for (auto& spline : splines) {
            for (int coord = 0; coord < 3; coord++) {
                const std::vector<double> spline_times = spline.GetTimes();

                const std::vector<int> position_nodes = spline.GetMutableNodes(EndEffectorSplines::Position, coord);
                for (auto& it: position_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 2.0;
                    spline.SetVars(EndEffectorSplines::Position, coord, it, vars);
                    REQUIRE(spline.ValueAt(EndEffectorSplines::Position, coord, spline_times.at(it)) == it);
                }

                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force, coord);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 2.0/FORCE_MULT;
                    spline.SetVars(EndEffectorSplines::Force, coord, it, vars);
                    REQUIRE(spline.ValueAt(EndEffectorSplines::Force, coord, spline_times.at(it)) == it);
                }
            }
        }
    }

    SECTION("Checking Values") {
        for (int coord = 0; coord < 3; coord++) {
            const std::vector<int> position_nodes = splines.at(0).GetMutableNodes(EndEffectorSplines::Position, coord);
            for (auto& it: position_nodes) {
                Eigen::Vector2d vars;
                vars << it, it - 1;
                splines.at(0).SetVars(EndEffectorSplines::Position, coord, it, vars);
            }
            REQUIRE(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0) == 0);
            if (coord != 2) {
                REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0.103448),
                             WithinAbs(1.0517, MARGIN));
                REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0.503448),
                             WithinAbs(4.62926, MARGIN));
            } else {
                REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0.162069),
                             WithinAbs(1.67752, MARGIN));
                REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Position, coord, 0.5),
                             WithinAbs(6, MARGIN));
            }
            const std::vector<int> position_nodes2 = splines.at(1).GetMutableNodes(EndEffectorSplines::Position, coord);
            for (auto& it: position_nodes2) {
                Eigen::Vector2d vars;
                vars << it, it - 1;
                splines.at(1).SetVars(EndEffectorSplines::Position, coord, it, vars);
            }
            REQUIRE(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0) == 0);
            if (coord != 2) {
                REQUIRE_THAT(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0.103448),
                             WithinAbs(0.0, MARGIN));
                REQUIRE_THAT(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0.25517),
                             WithinAbs(0.93156, MARGIN));
            } else {
                REQUIRE_THAT(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0.162069),
                             WithinAbs(0, MARGIN));
                REQUIRE_THAT(splines.at(1).ValueAt(EndEffectorSplines::Position, coord, 0.25517),
                             WithinAbs(2.2683, MARGIN));
            }
            const std::vector<int> force_nodes = splines.at(0).GetMutableNodes(EndEffectorSplines::Force, coord);
            for (auto& it: force_nodes) {
                Eigen::Vector2d vars;
                vars << it, (it - 1)/FORCE_MULT;
                splines.at(0).SetVars(EndEffectorSplines::Force, coord, it, vars);
            }
            REQUIRE(splines.at(0).ValueAt(EndEffectorSplines::Force, coord, 0) == 0);
            REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Force, coord, 0.103448), WithinAbs(0.0, MARGIN));
            REQUIRE_THAT(splines.at(0).ValueAt(EndEffectorSplines::Force, coord, 0.26666 + 0.0229885),
                         WithinAbs(3.27887, MARGIN));
        }
    }

    SECTION("Linearization/Coefficients") {
        for (auto& spline : splines) {
            for (int coord = 0; coord < 3; coord++) {
                const std::vector<double> spline_times = spline.GetTimes();
                const double total_time = spline.GetEndTime();

                const std::vector<int> position_nodes = spline.GetMutableNodes(EndEffectorSplines::Position, coord);
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

                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force, coord);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << it, 1.4/FORCE_MULT;
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
                REQUIRE_THAT(spline.GetEndTime(), WithinAbs(times.at(times.size() - 1) + (i + 1) * addt_time, MARGIN));
            }

            for (int coord = 0; coord < 3; coord++) {
                // Checking values after adding polys
                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force, coord);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << it - 1, .75/FORCE_MULT;
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
                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force, coord);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << 2 * it - 1, .5/FORCE_MULT;
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

                for (int i = 0; i < 20; i++) {
                    spline.RemovePoly(0.5 + i*0.1);
                }
            }
        }
    }

    SECTION("Derivatives") {
        for (auto& spline : splines) {
            // force
            for (int coord = 0; coord < 3; coord++) {
                const std::vector<int> force_nodes = spline.GetMutableNodes(EndEffectorSplines::Force, coord);
                for (auto& it: force_nodes) {
                    Eigen::Vector2d vars;
                    vars << 2 * it - 1, .5 / FORCE_MULT;
                    spline.SetVars(EndEffectorSplines::Force, coord, it, vars);
                }
            }

            time_v contact_times = spline.GetContactTimes();
            time_v contact_times2 = contact_times;
            EndEffectorSplines spline2 = spline;

            const double dt = std::sqrt(1e-16);
            constexpr double DERIV_MARGIN = 1e-4;

            // the last third part of the spline stuggles wrt to lift off

            double time = 0;
            while (time < spline.GetEndTime()) {
                for (int contact = 0; contact < contact_times.size(); contact++) {

                    for (int coord = 0; coord < 3; coord++) {
                        double val1 = spline.ValueAt(mpc::EndEffectorSplines::Force, coord, time);

                        contact_times2.at(contact).SetTime(contact_times2.at(contact).GetTime() + dt);
                        spline2.SetContactTimes(contact_times2);
                        double val2 = spline2.ValueAt(mpc::EndEffectorSplines::Force, coord, time);

                        double finite_diff = (val2 - val1)/(1*dt);

                        double partial = spline.ComputePartialWrtTime(mpc::EndEffectorSplines::Force, coord, time,
                                                                      contact);

                        REQUIRE_THAT(partial - finite_diff, WithinAbs(0, DERIV_MARGIN));

                        contact_times2.at(contact).SetTime(contact_times2.at(contact).GetTime() - dt);
                        spline2.SetContactTimes(contact_times2);
                    }
                }
                time += 0.01;
            }


            // Position
            for (int coord = 0; coord < 2; coord++) {
                const std::vector<int> pos_nodes = spline.GetMutableNodes(EndEffectorSplines::Position, coord);
                for (auto& it: pos_nodes) {
                    Eigen::Vector2d vars;
                    vars << 2 * it - 1, .5;
                    spline.SetVars(EndEffectorSplines::Position, coord, it, vars);
                }
            }

            spline2 = spline;

            // the last third part of the spline stuggles wrt to lift off

            time = 0;
            while (time < spline.GetEndTime()) {
                for (int contact = 0; contact < contact_times.size(); contact++) {

                    for (int coord = 0; coord < 2; coord++) {
                        double val1 = spline.ValueAt(mpc::EndEffectorSplines::Position, coord, time);

                        contact_times2.at(contact).SetTime(contact_times2.at(contact).GetTime() + dt);
                        spline2.SetContactTimes(contact_times2);
                        double val2 = spline2.ValueAt(mpc::EndEffectorSplines::Position, coord, time);

                        double finite_diff = (val2 - val1)/(1*dt);

                        double partial = spline.ComputePartialWrtTime(mpc::EndEffectorSplines::Position, coord, time,
                                                                      contact);

                        REQUIRE_THAT(partial - finite_diff, WithinAbs(0, DERIV_MARGIN));

                        contact_times2.at(contact).SetTime(contact_times2.at(contact).GetTime() - dt);
                        spline2.SetContactTimes(contact_times2);
                    }
                }
                time += 0.01;
            }
        }
    }
}