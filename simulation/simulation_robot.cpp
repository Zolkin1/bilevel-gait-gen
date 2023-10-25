//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "simulation_robot.h"

namespace simulator {
    SimulationRobot::SimulationRobot(const std::string& robot_xml_path, std::unique_ptr<Controller>& controller) :
        low_level_controller_(std::move(controller)) {
        robot_xml_path_ = robot_xml_path;
    }

    void SimulationRobot::SetInitialCondition(const Eigen::VectorXd& initial_config, const Eigen::VectorXd& initial_vel) {
        initial_config_ = initial_config;
        initial_vel_ = initial_vel;
    }

    std::string SimulationRobot::GetRobotXMLFile() const {
        return robot_xml_path_;
    }

    Eigen::VectorXd SimulationRobot::GetConfig() const {
        return initial_config_;
    }

    Eigen::VectorXd SimulationRobot::GetVelocities() const {
        return initial_vel_;
    }

    const Controller* SimulationRobot::GetController() const {
        return low_level_controller_.get();
    }

    void SimulationRobot::GetControlAction(const mjModel* model, mjtNum* cntrl) {
        std::vector<mjtNum> control = low_level_controller_->ComputeControlAction(model);
        for (int i = 0; i < model->nu; i++) {
            cntrl[i] = control.at(i);
        }
    }

} // simulator