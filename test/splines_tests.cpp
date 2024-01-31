//
// Created by zolkin on 1/30/24.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Core>
#include <iostream>

#include "end_effector_splines.h"

TEST_CASE("End Effector Splines", "[splines]") {
    using namespace mpc;
    std::vector<double> times;
    const int num_contacts = 10;
    const bool start_constant = true;
    const int num_force_polys = 3;

    for (int i = 0; i < num_contacts; i++) {
        times.emplace_back(0.2*i);
    }

    std::vector<mpc::EndEffectorSplines> splines;
    splines.emplace_back(num_contacts, times, start_constant, num_force_polys);
    splines.emplace_back(num_contacts, times, !start_constant, num_force_polys);


    SECTION("Setting Vars") {
        for (auto& spline : splines) {
            int temp = 0;
            for (int node = 0; node < spline.GetNumNodes(); node++) {
                Eigen::Vector2d vars;
                vars << temp, 2;

                if (spline.GetNodeType(EndEffectorSplines::Position, node) != NodeType::Empty) {
                    spline.SetVars(EndEffectorSplines::Position, node, vars);
                    temp++;
                }
                spline.SetVars(EndEffectorSplines::Force, node, vars);

//                std::cout << spline.ValueAt(EndEffectorSplines::Position, times.at(node)) << std::endl;
            }

            for (int node = 0; node < times.size()-1; node++) {
                REQUIRE(spline.ValueAt(EndEffectorSplines::Position, times.at(node)) == node);
                REQUIRE(spline.ValueAt(EndEffectorSplines::Force, times.at(node)) == node);
            }
        }
    }
}