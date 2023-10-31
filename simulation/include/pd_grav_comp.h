//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_PD_GRAV_COMP_H
#define BILEVEL_GAIT_GEN_PD_GRAV_COMP_H

#include "controller.h"

namespace simulator {
    class PDGravComp : public Controller {
    public:
        PDGravComp(double control_freq, std::string robot_urdf, const std::string& foot_type,
                   Eigen::VectorXd config_set_point, Eigen::VectorXd vel_set_point);

        /**
         * Computes the control action given the curren state and the set points.
         * When assigning values to the controllers, it is assumed that
         * 1. the yaml file aligns with the XML in terms of configuration and velocity
         * 2. the XML defines all the position actuators, then all the velocity acutators, then all the motors (feed forward)
         * @param model
         * @return a vector of control inputs to be given to the PD control
         */
        std::vector<mjtNum> ComputeControlAction(const mjModel* model, const mjData* data) override;

    private:
        /**
         * Computes the feed forward torque given the current state.
         * Assigns the feedforward torque to the private vector.
         */
        void ComputeFeedForwardValue(const mjData* data);

        /**
         * Assigns the configuration set point as the desired position
         * Note: assumes position actuators are first
         */
        void AssignPositionControl(std::vector<mjtNum>& control);

        /**
         * Assigns the velocity set point as the desired velocity
         * Note: assumes velocity actuators are second
         */
        void AssignVelocityControl(std::vector<mjtNum>& control);

        /**
         * Assigns the feedforward set point
         * Note: assumes feedforward acutators are last
         */
        void AssignFeedForward(std::vector<mjtNum>& control);


        Eigen::VectorXd config_set_point_;
        Eigen::VectorXd vel_set_point_;
        Eigen::VectorXd feedforward_;
    };
} // simulator

#endif //BILEVEL_GAIT_GEN_PD_GRAV_COMP_H
