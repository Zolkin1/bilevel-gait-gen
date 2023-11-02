//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CONTROLLER_H
#define BILEVEL_GAIT_GEN_CONTROLLER_H

#include <vector>

#include <pinocchio/multibody/data.hpp>
#include "pinocchio/parsers/urdf.hpp"
#include "mujoco.h"

namespace simulator {
    /**
     * Base class for low level controller for the robot.
     */
    class Controller {
    public:
        Controller(double control_freq, std::string robot_urdf, const std::string& foot_type);

        virtual void InitSolver(const mjData* data);

        /**
         * Specifies which joints we should expect to be in contact with the world.
         */
        void DefineContacts(const std::vector<std::string>& frames, const std::vector<int>& mujoco_bodies);

        double GetRate() const;

        int GetNumInputs() const;

        void PrintConfigNames() const;

        int GetNumContacts() const;

        /**
        * Interface with Mujoco to provide the current control action.
        */
        virtual std::vector<mjtNum> ComputeControlAction(const mjModel* model, const mjData* data) = 0;


        /**
         * Creates the mapping of the joints from mujoco to pinocchio.
         * Note: Assumes the joints have the same names in both the URDF and the XML.
         */
        void CreateJointMap(const mjModel* model);

    protected:
        Eigen::VectorXd ConvertMujocoConfigToPinocchio(const mjData* data) const;
        Eigen::VectorXd ConvertMujocoVelToPinocchio(const mjData* data) const;
        Eigen::VectorXd ConvertPinocchioJointToMujoco(const Eigen::VectorXd& joints);

        void UpdateContacts(const mjModel* model, const mjData* data);

        static constexpr int FLOATING_BASE_OFFSET = 7;
        static constexpr int FLOATING_VEL_OFFSET = 6;
        static constexpr int CONSTRAINT_PER_POINT_FOOT = 3;
        static constexpr int CONSTRAINT_PER_FLAT_FOOT = 6;

       int CONSTRAINT_PER_FOOT;

        double rate_;
        int num_inputs_;

        std::string robot_urdf_;

        pinocchio::Model pin_model_;
        std::unique_ptr<pinocchio::Data> pin_data_;

        std::vector<int> contact_frames_;       // frames potentially in contact
        std::vector<int> mujoco_bodies_;
        std::vector<bool> in_contact_;          // if each frame is in contact

        std::map<int, int> mujoco_to_pinocchio_joint_map_;
        std::vector<int> mujoco_joint_keys_;
    private:
    };
} // simulator


#endif //BILEVEL_GAIT_GEN_CONTROLLER_H
