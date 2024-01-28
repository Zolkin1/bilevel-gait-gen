//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "simulation_robot.h"

namespace simulator {
    SimulationRobot::SimulationRobot(const std::string& robot_xml_path, std::unique_ptr<controller::Controller>& controller) :
        low_level_controller_(std::move(controller)) {
        robot_xml_path_ = robot_xml_path;
        joint_map_created_ = false;
    }

    void SimulationRobot::SetInitialCondition(const Eigen::VectorXd& initial_config, const Eigen::VectorXd& initial_vel) {
        initial_config_ = initial_config;
        initial_vel_ = initial_vel;
    }

    void SimulationRobot::InitController() {
        low_level_controller_->InitSolver();
    }

    void SimulationRobot::InitController(const Eigen::VectorXd& full_body_state, const Eigen::VectorXd& mpc_state) {
        low_level_controller_->InitSolver(full_body_state, mpc_state);
    }

    void SimulationRobot::UpdateTargetConfig(const Eigen::VectorXd& q) {
        if (!joint_map_created_) {
            throw std::runtime_error("Can't assign target configuration until a joint map is created.");
        }
//        Eigen::VectorXd q_pin = ConvertMujocoVecConfigToPinocchio(q);
        low_level_controller_->UpdateTargetConfig(q);
    }

    void SimulationRobot::UpdateTargetVel(const Eigen::VectorXd& v) {
        if (!joint_map_created_) {
            throw std::runtime_error("Can't assign target velocity until a joint map is created.");
        }
//        Eigen::VectorXd v_pin = ConvertMujocoVecVelLikeToPinocchio(v);
        low_level_controller_->UpdateTargetVel(v);
    }

    void SimulationRobot::SetSimModel(const mjModel* model) {
        muj_model_ = model;
        num_inputs_ = muj_model_->nq - FLOATING_BASE_OFFSET;

        CreateJointMap();
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

    const controller::Controller* SimulationRobot::GetController() const {
        return low_level_controller_.get();
    }

    // TODO: make data const again
    void SimulationRobot::GetControlAction(const mjData* data, mjtNum* cntrl) {
        // Check that the dimensions allign
        if (muj_model_->nu != 3 * low_level_controller_->GetNumInputs()) {
            std::cerr << "Input mismatch! Mujoco is expecting " << muj_model_->nu/3 << " inputs while Pinocchio expects "
                      << low_level_controller_->GetNumInputs() << " inputs. Returning no control action." << std::endl;

            return;
        }

        Eigen::VectorXd q = ConvertMujocoConfigToPinocchio(data);
        Eigen::VectorXd v = ConvertMujocoVelToPinocchio(data);
        Eigen::VectorXd a = ConvertMujocoAccToPinocchio(data);

        // Update contact before querying controller
        UpdateContacts(data);

        // Compute control actions
        Eigen::VectorXd control = low_level_controller_->ComputeControlAction(q, v, a, contact_, data->time);

        // Convert the control back to mujoco
        std::vector<mjtNum> muj_control = ConvertControlToMujoco(control);

        for (int i = 0; i < muj_model_->nu; i++) {
            cntrl[i] = muj_control.at(i);
        }
    }

    // ---------------------- Joint Map ---------------------- //
    void SimulationRobot::CreateJointMap() {
        for (int i = 0; i < muj_model_->njnt; i++) {
            for (int j = 0; j < low_level_controller_->GetPinocchioModel().njoints; j++) {
                if (muj_model_->names + muj_model_->name_jntadr[i] ==
                low_level_controller_->GetPinocchioModel().names.at(j)) {
                    mujoco_to_pinocchio_joint_map_.insert(std::pair<int, int>(i,j));
                    mujoco_joint_keys_.push_back(i);
                }
            }
        }

        joint_map_created_ = true;
    }

    // ---------------------- Contacts ---------------------- //
    void SimulationRobot::DefineContacts(const std::vector<std::string>& frames, const std::vector<int>& mujoco_bodies) {
        if (frames.size() != mujoco_bodies.size()) {
            std::cerr << "Number of pinocchio contact frames and number of mujoco contact bodies do not match!"
                      << "Not assigning contacts in the controller." << std::endl;
            return;
        }

        pinocchio::Model pin_model = low_level_controller_->GetPinocchioModel();
        for (int i = 0; i < frames.size(); i++) {
            for (int j = 0; j < pin_model.frames.size(); j++) {
                // Check if the names given match a frame in pinocchio then take that index
                if (pin_model.frames.at(j).name == frames.at(i)) {
                    contact_.contact_frames_.push_back(j);
                    mujoco_bodies_.push_back(mujoco_bodies.at(i));
                    contact_.in_contact_.push_back(false);
                }
            }
        }

        // Check that we got everything
        if (contact_.contact_frames_.size() != frames.size()) {
            std::cerr << "Could not find " << frames.size() - contact_.contact_frames_.size()
                      << " frames. Check provided frames." << std::endl;
        }
    }

    int SimulationRobot::UpdateContacts(const mjData* data) {
        for (int i = 0; i < contact_.in_contact_.size(); i++) {
            contact_.in_contact_.at(i) = false;
        }

        int num_contacts = 0;

        for (int i = 0; i < data->ncon; i++) {
            for (int j = 0; j < mujoco_bodies_.size(); j++) {
                if (muj_model_->geom_bodyid[data->contact[i].geom2] == mujoco_bodies_.at(j)) {
                    contact_.in_contact_.at(j) = true;
                    num_contacts++;
                }
            }
        }

        return num_contacts;
    }

    int SimulationRobot::GetNumContacts() const {
        int num_contacts = 0;
        for (const bool contact : contact_.in_contact_) {
            if (contact) {
                num_contacts++;
            }
        }

        return num_contacts;
    }

    // ---------------------- Mujoco to Pinocchio Converters ---------------------- //
    Eigen::VectorXd SimulationRobot::ConvertMujocoConfigToPinocchio(const mjData* data) const {
        Eigen::VectorXd q = Eigen::VectorXd::Zero(low_level_controller_->GetPinocchioModel().nq);

        // Floating base config
        // floating base position
        for (int i = 0; i < 3; i++) {
            q(i) = data->qpos[i];
        }

        // floating base quaternion, note pinocchio uses (x,y,z,w_) and mujoco uses (w_,x,y,z)
        q(6) = data->qpos[3];
        q(3) = data->qpos[4];
        q(4) = data->qpos[5];
        q(5) = data->qpos[6];


        // Joints
        for (int i = 0; i < mujoco_joint_keys_.size(); i++) {
            q(mujoco_to_pinocchio_joint_map_.at(mujoco_joint_keys_.at(i)) - 2 + FLOATING_BASE_OFFSET) = data->qpos[i + FLOATING_BASE_OFFSET];
        }


        return q;
    }

    Eigen::VectorXd SimulationRobot::ConvertMujocoVelToPinocchio(const mjData* data) const {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(low_level_controller_->GetPinocchioModel().nv);

        // Floating base velocities
        for (int i = 0; i < FLOATING_VEL_OFFSET; i++) {
            v(i) = data->qvel[i];
        }

        // Joint velocities
        for (int i = 0; i < mujoco_joint_keys_.size(); i++) {
            v(mujoco_to_pinocchio_joint_map_.at(mujoco_joint_keys_.at(i)) - 2 + FLOATING_VEL_OFFSET) = data->qvel[i + FLOATING_VEL_OFFSET];
        }

        return v;
    }

    Eigen::VectorXd SimulationRobot::ConvertMujocoAccToPinocchio(const mjData* data) const {
        Eigen::VectorXd a = Eigen::VectorXd::Zero(low_level_controller_->GetPinocchioModel().nv);

        // Floating base velocities
        for (int i = 0; i < FLOATING_VEL_OFFSET; i++) {
            a(i) = data->qacc[i];
        }

        // Joint velocities
        for (int i = 0; i < mujoco_joint_keys_.size(); i++) {
            a(mujoco_to_pinocchio_joint_map_.at(mujoco_joint_keys_.at(i)) - 2 + FLOATING_VEL_OFFSET) =
                    data->qacc[i + FLOATING_VEL_OFFSET];
        }

        return a;
    }

    Eigen::VectorXd SimulationRobot::ConvertMujocoVecConfigToPinocchio(const Eigen::VectorXd& q) const {
        Eigen::VectorXd q_out = Eigen::VectorXd::Zero(low_level_controller_->GetPinocchioModel().nq);

        // Floating base config
        // floating base position
        for (int i = 0; i < 3; i++) {
            q_out(i) = q(i);
        }

        // floating base quaternion, note pinocchio uses (x,y,z,w_) and mujoco uses (w_,x,y,z)
        q_out(6) = q(3);
        q_out(3) = q(4);
        q_out(4) = q(5);
        q_out(5) = q(6);


        // Joints
        for (int i = 0; i < mujoco_joint_keys_.size(); i++) {
            q_out(mujoco_to_pinocchio_joint_map_.at(mujoco_joint_keys_.at(i)) - 2 + FLOATING_BASE_OFFSET) =
                    q(i + FLOATING_BASE_OFFSET);
        }

        return q_out;
    }

    Eigen::VectorXd SimulationRobot::ConvertMujocoVecVelLikeToPinocchio(const Eigen::VectorXd& v) const {
        Eigen::VectorXd v_out = Eigen::VectorXd::Zero(low_level_controller_->GetPinocchioModel().nv);

        // Floating base velocities
        for (int i = 0; i < FLOATING_VEL_OFFSET; i++) {
            v_out(i) = v(i);
        }

        // Joint velocities
        for (int i = 0; i < mujoco_joint_keys_.size(); i++) {
            v_out(mujoco_to_pinocchio_joint_map_.at(mujoco_joint_keys_.at(i)) - 2 + FLOATING_VEL_OFFSET) =
                    v(i + FLOATING_VEL_OFFSET);
        }

        return v_out;
    }


    // ---------------------- Pinocchio to Mujoco Converters ---------------------- //
    Eigen::VectorXd SimulationRobot::ConvertPinocchioJointToMujoco(const Eigen::VectorXd& joints) const {
        Eigen::VectorXd mujoco_joints(joints.size());

        for (int i = 0; i < mujoco_joint_keys_.size(); i++) {
            mujoco_joints(i) = joints(mujoco_to_pinocchio_joint_map_.at(mujoco_joint_keys_.at(i)) - 2);
        }

        return mujoco_joints;
    }

    Eigen::VectorXd SimulationRobot::ConvertPinocchioVelToMujoco(const Eigen::VectorXd& v) const {
        Eigen::VectorXd mujoco_vel(v.size());
        // floating base
        for (int i = 0; i < FLOATING_VEL_OFFSET; i++) {
            mujoco_vel(i) = v(i);
        }

        // joints
        mujoco_vel.tail(num_inputs_) = ConvertPinocchioJointToMujoco(v.tail(num_inputs_));

        return mujoco_vel;
    }

    Eigen::VectorXd SimulationRobot::ConvertPinocchioConfigToMujoco(const Eigen::VectorXd& q) const{
        assert(q.size() == muj_model_->nq);
        Eigen::VectorXd muj_q(q.size());
        muj_q.tail(muj_model_->nv - FLOATING_VEL_OFFSET) =
                ConvertPinocchioJointToMujoco(q.tail(muj_model_->nv - FLOATING_VEL_OFFSET));
        muj_q.head(3) = q.head(3);

        // floating base quaternion, note pinocchio uses (x,y,z,w_) and mujoco uses (w_,x,y,z)
        muj_q(3) = q(6);
        muj_q(5) = q(4);
        muj_q(6) = q(5);
        muj_q(4) = q(3);
        return muj_q;
    }

    std::vector<mjtNum> SimulationRobot::ConvertControlToMujoco(const Eigen::VectorXd& control) const {
        std::vector<mjtNum> muj_control(3*control.size());
        for (int i = 0; i < 3; i++) {
            Eigen::VectorXd control_part = ConvertPinocchioJointToMujoco(control.segment(i*num_inputs_, num_inputs_));
            for (int j = 0; j < num_inputs_; j++) {
                muj_control.at(j + i*num_inputs_) = control_part(j);
            }
        }

        return muj_control;
    }

    std::vector<std::vector<Eigen::Vector3d>> SimulationRobot::GetTrajViz() {
        return low_level_controller_->GetTrajViz();
    }

} // simulator