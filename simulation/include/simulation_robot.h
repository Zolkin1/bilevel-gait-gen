//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H
#define BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H

#include <string>
#include <Eigen/Core>
#include <memory>

#include "mujoco.h"
#include <pinocchio/multibody/data.hpp>

#include "controller.h"

/**
 * SimulationRobot class.
 * Defines a robot to be used by the simulator. Interfaces with the controller.
 */
namespace simulator {

    // TODO: This class should really handle ALL of the mujoco interface so this interface can be swapped for hardware.
    class SimulationRobot {
    public:
        SimulationRobot(const std::string& robot_xml_path, std::unique_ptr<controller::Controller>& controller);

        void SetInitialCondition(const Eigen::VectorXd& initial_condition, const Eigen::VectorXd& initial_vel);

        void InitController();

        void UpdateTargetConfig(const Eigen::VectorXd& q);

        void UpdateTargetVel(const Eigen::VectorXd& v);

        void SetSimModel(const mjModel* model);

        [[nodiscard]] std::string GetRobotXMLFile() const;

        [[nodiscard]] Eigen::VectorXd GetInitConfig() const;

        [[nodiscard]] Eigen::VectorXd GetInitVelocities() const;

        /**
         * Get the low level controller
         * @return a const pointer to the controller
         */
        [[nodiscard]] const controller::Controller* GetController() const;

        /**
         * Interface with Mujoco to provide the current control action.
         * Puts the controls into the control vector provided.
         */
        void GetControlAction(const mjData* data, mjtNum* control);

        /**
         * Creates the mapping of the joints from mujoco to pinocchio.
         * Note: Assumes the joints have the same names in both the URDF and the XML.
         */
        void CreateJointMap();

        /**
         * Specifies which joints we should expect to be in contact with the world.
         */
        void DefineContacts(const std::vector<std::string>& frames, const std::vector<int>& mujoco_bodies);
        int UpdateContacts(const mjData* data);
        int GetNumContacts() const;

        Eigen::VectorXd ConvertMujocoConfigToPinocchio(const mjData* data) const;
        Eigen::VectorXd ConvertMujocoVelToPinocchio(const mjData* data) const;
        Eigen::VectorXd ConvertMujocoAccToPinocchio(const mjData* data) const;

        Eigen::VectorXd ConvertMujocoVecConfigToPinocchio(const Eigen::VectorXd& q) const;
        Eigen::VectorXd ConvertMujocoVecVelLikeToPinocchio(const Eigen::VectorXd& v) const;

        Eigen::VectorXd ConvertPinocchioJointToMujoco(const Eigen::VectorXd& joints);
        Eigen::VectorXd ConvertPinocchioVelToMujoco(const Eigen::VectorXd& v);
        std::vector<mjtNum> ConvertControlToMujoco(const Eigen::VectorXd& control);

    private:
        std::string robot_xml_path_;
        std::unique_ptr<controller::Controller> low_level_controller_;

        Eigen::VectorXd initial_config_;
        Eigen::VectorXd initial_vel_;

        bool joint_map_created_;

        Eigen::VectorXd target_config_;
        Eigen::VectorXd target_vel_;

        static constexpr int FLOATING_BASE_OFFSET = 7;
        static constexpr int FLOATING_VEL_OFFSET = 6;

        // Mappings between mujoco and controller interface
        std::map<int, int> mujoco_to_pinocchio_joint_map_;
        std::vector<int> mujoco_joint_keys_;

        int num_inputs_;

        // Internal copy of the mujoco model
        const mjModel* muj_model_;

        // Contact information
//        std::vector<int> contact_frames_;       // frames potentially in contact
        std::vector<int> mujoco_bodies_;
//        std::vector<bool> in_contact_;          // if each frame is in contact

        controller::Contact contact_;

    };
} // simulator


#endif //BILEVEL_GAIT_GEN_SIMULATION_ROBOT_H
