//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_CONTROLLER_H
#define BILEVEL_GAIT_GEN_CONTROLLER_H

#include <vector>

#include <pinocchio/multibody/data.hpp>
#include "pinocchio/parsers/urdf.hpp"

namespace controller {

    struct Contact {
        std::vector<bool> in_contact_;          // if each frame is in contact
        std::vector<int> contact_frames_;       // frames potentially in contact

        Contact();

        Contact(int num_contacts);

        int GetNumContacts() const;
    };

    /**
     * Base class for low level controller for the robot.
     */
    class Controller {
    public:
        Controller(double control_freq, std::string robot_urdf, const std::string& foot_type);

        virtual void InitSolver();

        double GetRate() const;

        int GetNumInputs() const;

        void PrintConfigNames() const;

        const pinocchio::Model& GetPinocchioModel() const;

        virtual Eigen::VectorXd ComputeControlAction(const Eigen::VectorXd& q,
                                                     const Eigen::VectorXd& v,
                                                     const Eigen::VectorXd& a,
                                                     const Contact& contact,
                                                     double time) = 0;

        void UpdateTargetConfig(const Eigen::VectorXd& q);
        void UpdateTargetVel(const Eigen::VectorXd& v);
        void UpdateTargetAcc(const Eigen::VectorXd& a);

        virtual void InitSolver(const Eigen::VectorXd& state);

        virtual void UpdateDesiredContacts(const Contact& contact);

        virtual std::vector<std::vector<Eigen::Vector3d>> GetTrajViz();

    protected:
        /**
         * Assigns the configuration set point as the desired position
         * Note: assumes position actuators are first
         */
        void AssignPositionControl(Eigen::VectorXd& control);

        /**
         * Assigns the velocity set point as the desired velocity
         * Note: assumes velocity actuators are second
         */
        void AssignVelocityControl(Eigen::VectorXd& control);

        static constexpr int FLOATING_BASE_OFFSET = 7;
        static constexpr int FLOATING_VEL_OFFSET = 6;
        static constexpr int CONSTRAINT_PER_POINT_FOOT = 3;
        static constexpr int CONSTRAINT_PER_FLAT_FOOT = 6;

        // Number of position variables
        static constexpr int POS_VARS = 3;

        int CONSTRAINT_PER_FOOT;

        double rate_;
        int num_inputs_;

        std::string robot_urdf_;

        pinocchio::Model pin_model_;
        std::unique_ptr<pinocchio::Data> pin_data_;

        Eigen::VectorXd config_target_;
        Eigen::VectorXd vel_target_;
        Eigen::VectorXd acc_target_;
    private:
    };
} // controller


#endif //BILEVEL_GAIT_GEN_CONTROLLER_H
