//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_PD_GRAV_COMP_H
#define BILEVEL_GAIT_GEN_PD_GRAV_COMP_H

#include "controller.h"

namespace controller {
    class PDGravComp : public Controller {
    public:
        PDGravComp(double control_freq, std::string robot_urdf, const std::string& foot_type,
                   const Eigen::VectorXd& config_set_point, const Eigen::VectorXd& vel_set_point);

        /**
         * Computes the control action given the curren state and the set points.
         * When assigning values to the controllers, it is assumed that
         * 1. the yaml file aligns with the XML in terms of configuration and velocity
         * 2. the XML defines all the position actuators, then all the velocity acutators, then all the motors (feed forward)
         * @param model
         * @return a vector of control inputs to be given to the PD control
         */
        Eigen::VectorXd ComputeControlAction(const Eigen::VectorXd& q,
                                                 const Eigen::VectorXd& v,
                                                 const Eigen::VectorXd& a,
                                                 const Contact& contact,
                                                 double time) override;

    private:
        /**
         * Computes the feed forward torque given the current state.
         * Assigns the feedforward torque to the private vector.
         */
        void ComputeFeedForwardValue(const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Contact& contact);

        /**
         * Assigns the feedforward set point
         * Note: assumes feedforward acutators are last
         */
        void AssignFeedForward(Eigen::VectorXd& control);

    };
} // controller

#endif //BILEVEL_GAIT_GEN_PD_GRAV_COMP_H
