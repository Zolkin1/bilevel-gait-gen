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

    void SimulationRobot::InitController(const mjModel* model, const mjData* data) {
        low_level_controller_->InitSolver(model, data);
    }

    std::string SimulationRobot::GetRobotXMLFile() const {
        return robot_xml_path_;
    }

    Eigen::VectorXd SimulationRobot::GetInitConfig() const {
        return initial_config_;
    }

    Eigen::VectorXd SimulationRobot::GetInitVelocities() const {
        return initial_vel_;
    }

    const Controller* SimulationRobot::GetController() const {
        return low_level_controller_.get();
    }

    void SimulationRobot::GetControlAction(const mjModel* model, const mjData* data, mjtNum* cntrl) {
        // Check that the dimensions allign
        if (model->nu != 3 * low_level_controller_->GetNumInputs()) {
            std::cerr << "Input mismatch! Mujoco is expecting " << model->nu/3 << " inputs while Pinocchio expects "
                      << low_level_controller_->GetNumInputs() << " inputs. Returning no control action." << std::endl;

            return;
        }

        // Compute control actions
        std::vector<mjtNum> control = low_level_controller_->ComputeControlAction(model, data);
        for (int i = 0; i < model->nu; i++) {
            cntrl[i] = control.at(i);
        }
    }

    void SimulationRobot::CreateJointMap(const mjModel* model) {
        low_level_controller_->CreateJointMap(model);
    }

} // simulator