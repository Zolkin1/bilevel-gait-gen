//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_PD_GRAV_COMP_H
#define BILEVEL_GAIT_GEN_PD_GRAV_COMP_H

#include "controller.h"

namespace simulator {
    class PDGravComp : public Controller {
    public:
        PDGravComp(double control_freq, Eigen::VectorXd set_point);

        std::vector<mjtNum> ComputeControlAction(const mjModel *model) override;

    private:
        Eigen::VectorXd set_point_;
    };
} // simulator

#endif //BILEVEL_GAIT_GEN_PD_GRAV_COMP_H
