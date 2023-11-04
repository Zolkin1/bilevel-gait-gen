//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H
#define BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H

#include <string>
#include <Eigen/Core>
#include <memory>

#include "mujoco.h"

#include "controller.h"

/**
 * SimulationRobot class.
 * Defines a robot to be used by the simulator. Holds all the information needed by the simulator.
 */
namespace simulator {

    // TODO: This class should really handle ALL of the mujoco interface so this interface can be swapped for hardware.
    class SimulationRobot {
    public:
        SimulationRobot(const std::string& robot_xml_path, std::unique_ptr<Controller>& controller);

        void SetInitialCondition(const Eigen::VectorXd& initial_condition, const Eigen::VectorXd& initial_vel);

        void InitController(const mjModel* model, const mjData* data);

        void UpdateTargetConfig(Eigen::VectorXd q);

        void UpdateTargetVel(Eigen::VectorXd v);

        [[nodiscard]] std::string GetRobotXMLFile() const;

        [[nodiscard]] Eigen::VectorXd GetInitConfig() const;

        [[nodiscard]] Eigen::VectorXd GetInitVelocities() const;

        /**
         * Get the low level controller
         * @return a const pointer to the controller
         */
        [[nodiscard]] const Controller* GetController() const;

        /**
         * Interface with Mujoco to provide the current control action.
         * Puts the controls into the control vector provided.
         */
        void GetControlAction(const mjModel* model, const mjData* data, mjtNum* control);

        void CreateJointMap(const mjModel* model);
    private:
        std::string robot_xml_path_;
        std::unique_ptr<Controller> low_level_controller_;

        Eigen::VectorXd initial_config_;
        Eigen::VectorXd initial_vel_;

        bool joint_map_created_;

        Eigen::VectorXd target_config_;
        Eigen::VectorXd target_vel_;
    };
} // simulator


#endif //BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H
