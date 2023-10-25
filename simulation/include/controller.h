//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CONTROLLER_H
#define BILEVEL_GAIT_GEN_CONTROLLER_H

#include <vector>

#include "mujoco.h"

namespace simulator {
    /**
     * Base class for low level controller for the robot.
     */
    class Controller {
    public:
        Controller(double control_freq);

        double GetRate() const;

        /**
        * Interface with Mujoco to provide the current control action.
        */
        virtual std::vector<mjtNum> ComputeControlAction(const mjModel* model) = 0;

    private:
        double rate_;
    };
} // simulator


#endif //BILEVEL_GAIT_GEN_CONTROLLER_H
