//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <Eigen/Core>

#include "mujoco.h"

#include "pd_grav_comp.h"


namespace simulator {
    PDGravComp::PDGravComp(double control_freq, Eigen::VectorXd set_point):Controller(control_freq),
    set_point_(set_point) {}

    std::vector<mjtNum> PDGravComp::ComputeControlAction(const mjModel* model) {
        std::vector<mjtNum> control;
        for (int i = 0; i < model->nu; i++) {
            if (i < 11) {       // TODO: Make this not hard coded
                control.push_back(set_point_(i + 7));
            } else {
                control.push_back(0);
            }
        }
        return control;
    }
}