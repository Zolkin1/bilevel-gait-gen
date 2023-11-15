//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/centroidal-derivatives.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/spatial/explog-quaternion.hpp"

#include "centroidal_model.h"
#include "euler_integrator.h"


namespace mpc {
    using matrix6x_t = Eigen::Matrix<double, 6, Eigen::Dynamic>;

    CentroidalModel::CentroidalModel(const std::string &robot_urdf, const std::vector<std::string>& frames,
                                     int discretization_steps, double dt) :
            GRAVITY(0., 0., -9.81) {

        integrator_ = std::make_unique<EulerIntegrator>(dt);

        // create the pinocchio model - always a free flyer
        pinocchio::urdf::buildModel(robot_urdf, pinocchio::JointModelFreeFlyer(), pin_model_, false);

        // create the pinocchio data
        pin_data_ = std::make_unique<pinocchio::Data>(pin_model_);

        num_joints_ = pin_model_.nq - FLOATING_BASE_OFFSET;
        num_total_states_ = MOMENTUM_OFFSET + FLOATING_VEL_OFFSET + num_joints_;

        robot_mass_ = pin_data_->mass[0];

        num_ee_ = frames.size();
        frames_= frames;

        CreateFrameMap(frames);

        if (discretization_steps != 1) {
            throw std::runtime_error("Only discretization step of 1 is currently supported.");
        }
        discretization_steps_ = discretization_steps;
    }

    void CentroidalModel::GetLinearDiscreteDynamics(const vector_t& state, const vector_t& ref_state, const Inputs& input,
                                                             double time, matrix_t& A, matrix_t& B,
                                                             vector_t& C) {

        // Get the pinocchio configuration
        vector_t q_pin = ConvertMPCStateToPinocchioState(state);

        // ------------------------------------------- //
        // --------------- Calculate A --------------- //
        // ------------------------------------------- //
        A = matrix_t::Zero(num_total_states_, num_total_states_);

        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q_pin);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q_pin);

        Eigen::Matrix3Xd cross_sum = Eigen::Matrix3Xd::Zero(3, pin_model_.nv);
        for (int i = 0; i < num_ee_; i++) {
            // Get the jacobian of the FK
            matrix_t J = GetFKJacobianForEndEffector(q_pin, frames_.at(i), false);
            Eigen::Vector3d force = input.GetForce(i, time);
            for (int j = 0; j < pin_model_.nv; j++) {
                // Cross product with the corresponding input force, add these for all the end effectors
                cross_sum.col(j) += static_cast<Eigen::Vector3d>(-J.col(j)).cross(force);
            }
        }
        A.block(3, MOMENTUM_OFFSET, cross_sum.rows(), cross_sum.cols()) = cross_sum;

        // Deal with centroidal dynamics
        matrix6x_t CMM = pinocchio::computeCentroidalMap(pin_model_, *pin_data_, q_pin);

        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
        matrix_t CMMj = CMM.rightCols(num_joints_);
        matrix_t CMMbinv = CMMb.inverse();  // TODO: maybe calc better
        vector_t v_pin = ComputePinocchioVelocities(state, input.GetVels(time), CMMbinv, CMMj);

        vector_t a_pin = vector_t::Zero(pin_model_.nv);

        matrix6x_t dhdq = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
        matrix6x_t dhdotdq = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
        matrix6x_t dhdotdv = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);

        // Get dhdq
        // TODO: I only need this for dhdq so I can just steal that one line for efficiency
        pinocchio::computeCentroidalDynamicsDerivatives(pin_model_, *pin_data_, q_pin, v_pin, a_pin,
                                                        dhdq, dhdotdq, dhdotdv, CMM);

        // TODO: Check dhdq
        A.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) << CMMbinv, CMMbinv * dhdq; // TODO: Check negative sign


        // TODO: For now I'm just assuming the CMM has no partial wrt q. Will need to come back and change this.
        // I need the partial w.r.t to each configuration variable
        // Note: (A^-1)' = (A^-1)*A'*(A^-1)

        // ------------------------------------------- //
        // --------------- Calculate B --------------- //
        // ------------------------------------------- //
        int num_inputs = input.GetNumInputs();

        // Inputs: [f_{c1,p1}^{x}(1), ... , f_{c1, p1}^{x}(4), f_{c1,p2}^{x}(1), ... , f_{c1, p2}^{x}(4), ... ]
        //          = [f_{c1,p1}^{x}, f_{c1,p2}^{x}, f_{c1,p3}^{x}, ... ]
        //          = [f_{c1}^{x}, f_{c1}^{y}, f_{c1}^{z}, f_{c2}^{x}, f_{c2}^{y}, f_{c2}^{z}, ... ]
        //     Given a specific time, we only need to look at one pj (polynomial j) and thus
        //     one set of columns for each spline.

        B = matrix_t::Zero(num_total_states_, num_inputs);
        // --- Linear momentum --- //
        // Linear in each force, and each force is linear in its coefficients at each point in time.
        int idx = input.GetForcePolyIdx(time);
        int vars_per_spline = input.GetForces().at(0).at(0).GetTotalPolyVars();
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < 3; coord++) {
                Eigen::Matrix<double, 3, 4> vars_lin = input.GetForcePolyVarsLin(ee, time);
                B.block<1, Spline::POLY_ORDER>(coord, ee*3*vars_per_spline + coord*vars_per_spline + idx*Spline::POLY_ORDER)
                        = vars_lin.row(coord);
            }
        }

        Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();
        // --- Angular momentum --- //
        for (int ee = 0; ee < num_ee_; ee++) {
            Eigen::Vector3d ee_pos_wrt_com = GetEndEffectorLocationCOMFrame(state, frames_.at(ee));
            for (int coord = 0; coord < 3; coord++) {
                Eigen::Matrix<double, 3, 4> vars_lin = input.GetForcePolyVarsLin(ee, time);
                for (int poly = 0; poly < Spline::POLY_ORDER; poly++) {
                    B.block<3, 1>(3, ee*3*vars_per_spline + coord*vars_per_spline + idx*Spline::POLY_ORDER + poly) =
                            ee_pos_wrt_com.cross(static_cast<Eigen::Vector3d>(Id.col(coord))) * vars_lin(coord, poly);
                }
            }
        }

        // --- Base velocity --- //
        B.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) <<
                matrix_t::Zero(FLOATING_VEL_OFFSET, num_inputs - num_joints_),
                -CMMbinv*CMMj*matrix_t::Identity(num_joints_, num_joints_);

        // --- Joint velocity --- //
        B.bottomRows(num_joints_) << matrix_t::Zero(num_joints_, num_inputs - num_joints_),
                matrix_t::Identity(num_joints_, num_joints_);


        // ------------------------------------------- //
        // ---------------- Calculate C -------------- //
        // ------------------------------------------- //
        C = CalcDynamics(state, input, time);

        // Discretize with the integrator
        A = integrator_->CalcDerivWrtStateSingleStep(state, A);
        B = integrator_->CalcDerivWrtInputSingleStep(state, B);

    }

    void CentroidalModel::GetFKLinearization(const vector_t& state, const Inputs& input, int end_effector,
                                                 matrix_t& A, vector_t& C) {
        vector_t q_pin = ConvertMPCStateToPinocchioState(state);
        A = matrix_t::Zero(3, pin_model_.nv);
        C = vector_t::Zero(3);

        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q_pin);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q_pin);

        matrix6x_t J = matrix6x_t::Zero(6, pin_model_.nv);
        pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_map_.at(frames_.at(end_effector)),
                                    pinocchio::LOCAL_WORLD_ALIGNED, J);

        A = J.topRows<3>();

        C = GetEndEffectorLocationCOMFrame(state, frames_.at(end_effector));
    }

    matrix_t CentroidalModel::GetFKJacobianForEndEffector(const vector_t& q, const std::string& frame, bool compute_jac) {
        if (compute_jac) {
            pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
            pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);
        }

        matrix6x_t J = matrix_t::Zero(6,pin_model_.nv);
        pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_map_.at(frame), pinocchio::LOCAL_WORLD_ALIGNED, J);

        // Only return the position coordinates
        return J.topRows<3>();
    }

    Eigen::Vector3d CentroidalModel::GetEndEffectorLocationCOMFrame(const vector_t& state, const std::string& frame) const {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state));
        return pin_data_->oMf.at(frame_map_.at(frame)).translation();
    }

    // Note: The MPC state is using a quaternion representation
    int CentroidalModel::GetPinocchioNumConfig() const {
        return pin_model_.nq;
    }

    int CentroidalModel::GetNumJoints() const {
        return num_joints_;
    }

    int CentroidalModel::GetNumEndEffectors() const {
        return num_ee_;
    }

    vector_t CentroidalModel::ComputePinocchioVelocities(const vector_t& state, const vector_t& joint_vels,
                                                         const matrix_t& CMMbinv, const matrix_t& CMMj) const {
        vector_t v_pin(FLOATING_VEL_OFFSET + num_joints_);

        v_pin.head<FLOATING_VEL_OFFSET>() = CMMbinv * (state.head<MOMENTUM_OFFSET>() - CMMj*joint_vels);

        v_pin.tail(num_joints_) = joint_vels;

        return v_pin;
    }

    vector_t CentroidalModel::ConvertMPCStateToPinocchioState(const vector_t &state) const {
        vector_t q(pin_model_.nq);

        q.head(POS_VARS) = state.segment<POS_VARS>(MOMENTUM_OFFSET);
        q.segment<4>(POS_VARS) = state.segment<4>(MOMENTUM_OFFSET + POS_VARS);
        q.tail(num_joints_) = state.tail(num_joints_);

        return q;
    }

    vector_t CentroidalModel::ConvertPinocchioStateToMPCState(const Eigen::Vector<double, MOMENTUM_OFFSET>& momentum,
                                                              const vector_t& state) const {
        Eigen::Vector3d rot = ConvertQuaternionToZYXRot(state.segment<4>(POS_VARS));
        vector_t x(num_total_states_);
        x.head<MOMENTUM_OFFSET>() = momentum;
        x.segment<POS_VARS>(MOMENTUM_OFFSET) = state.head<POS_VARS>();
        x.segment<3>(MOMENTUM_OFFSET + POS_VARS) = rot;
        x.tail(num_joints_) = state.tail(num_joints_);

        return x;
    }

    // converts to a x,y,z,w quaternion
    Eigen::Vector4d CentroidalModel::ConvertZYXRotToQuaternion(const Eigen::Vector3d& zyx_rot) {
        Eigen::Quaterniond quat = static_cast<Eigen::Quaterniond>(Eigen::AngleAxisd(zyx_rot(0), Eigen::Matrix<double, 3, 1>::UnitZ()) *
                                           Eigen::AngleAxisd(zyx_rot(1), Eigen::Matrix<double, 3, 1>::UnitY()) *
                                           Eigen::AngleAxisd(zyx_rot(2), Eigen::Matrix<double, 3, 1>::UnitX()));
        Eigen::Vector4d quat_vec;
        quat_vec(0) = quat.x();
        quat_vec(1) = quat.y();
        quat_vec(2) = quat.z();
        quat_vec(3) = quat.w();
        return quat_vec;
    }

    Eigen::Vector3d CentroidalModel::ConvertQuaternionToZYXRot(const Eigen::Vector4d& quat) {
        Eigen::Vector3d rot;
        double q_w = quat(3);
        double q_x = quat(0);
        double q_y = quat(1);
        double q_z = quat(2);
        // x
        rot(2) = std::atan2(2*(q_w*q_x + q_y*q_z), 1 - 2*(q_x*q_x + q_y*q_y));
        // y
        rot(1) = -M_PI/2 + 2*std::atan2(std::sqrt(1 + 2*(q_w*q_y - q_x*q_z)), std::sqrt(1 - 2*(q_w*q_y - q_x*q_z)));
        // z
        rot(0) = std::atan2(2*(q_w*q_z + q_x*q_y), 1 - 2*(q_y*q_y + q_z*q_z));

        return rot;
    }

    vector_t CentroidalModel::ConvertManifoldStateToAlgebraState(const vector_t& state, const vector_t& ref_state) {
        vector_t alg_state(state.size()-1);
        alg_state.head<MOMENTUM_OFFSET>() = state.head<MOMENTUM_OFFSET>();
        alg_state.segment<POS_VARS>(MOMENTUM_OFFSET) = state.segment<POS_VARS>(MOMENTUM_OFFSET);
        Eigen::Quaterniond quat(static_cast<Eigen::Vector4d>(state.segment<4>(POS_VARS + MOMENTUM_OFFSET)));
        Eigen::Quaterniond quat_ref(static_cast<Eigen::Vector4d>(ref_state.segment<4>(POS_VARS + MOMENTUM_OFFSET)));
        alg_state.segment<3>(POS_VARS + MOMENTUM_OFFSET) = pinocchio::quaternion::log3(quat_ref.inverse()*quat);
        alg_state.tail(state.size() - (MOMENTUM_OFFSET + FLOATING_BASE_OFFSET)) =
                state.tail(state.size() - (MOMENTUM_OFFSET + FLOATING_BASE_OFFSET));
        return alg_state;
    }

    vector_t CentroidalModel::ConvertAlgebraStateToManifoldState(const mpc::vector_t &state, const vector_t& ref_state) {
        vector_t man_state(state.size()+1);
        Eigen::Quaterniond quat;
        pinocchio::quaternion::exp3(state.segment<3>(MOMENTUM_OFFSET+POS_VARS), quat);
        Eigen::Quaterniond quat_ref(static_cast<Eigen::Vector4d>(ref_state.segment<4>(POS_VARS + MOMENTUM_OFFSET)));
        quat = quat_ref*quat;
        man_state(0+MOMENTUM_OFFSET+POS_VARS) = quat.x();
        man_state(1+MOMENTUM_OFFSET+POS_VARS) = quat.y();
        man_state(2+MOMENTUM_OFFSET+POS_VARS) = quat.z();
        man_state(3+MOMENTUM_OFFSET+POS_VARS) = quat.w();

        man_state.head<MOMENTUM_OFFSET>() = state.head<MOMENTUM_OFFSET>();
        man_state.segment<POS_VARS>(MOMENTUM_OFFSET) = state.segment<POS_VARS>(MOMENTUM_OFFSET);
        man_state.tail(state.size()-MOMENTUM_OFFSET-FLOATING_VEL_OFFSET) =
                state.tail(state.size()-MOMENTUM_OFFSET-FLOATING_VEL_OFFSET);

        return man_state;
    }

    void CentroidalModel::CreateFrameMap(const std::vector<std::string>& frames) {
        for (int i = 0; i < num_ee_; i++) {
            for (int j = 0; j < pin_model_.frames.size(); j++)
            if (pin_model_.frames.at(j).name == frames.at(i)) {
                frame_map_.insert(std::pair<std::string, int>(frames.at(i), pin_model_.getFrameId(frames.at(i))));
            }
        }
    }

    vector_t CentroidalModel::CalcDynamics(const vector_t& state, const Inputs& input, double time) const {
        vector_t xdot = vector_t::Zero(num_total_states_);

        xdot.topRows<3>() = robot_mass_*GRAVITY;

        for (int i = 0; i < num_ee_; i++) {
            xdot.topRows<3>() += input.GetForce(i, time);

            // Note that the input ee position can deviate from the state ee position.
            // So I use the state ee position.
            xdot.middleRows<3>(3) += GetEndEffectorLocationCOMFrame(state, frames_.at(i)).cross(input.GetForce(i, time));
        }

        matrix_t CMM = pinocchio::computeCentroidalMap(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state));
        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
        matrix_t CMMj = CMM.rightCols(num_joints_);
        xdot.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) =
                CMMb.householderQr().solve(state.head(MOMENTUM_OFFSET) - CMMj*input.GetVels(time));

        xdot.bottomRows(num_joints_) = input.GetVels(time);

        return xdot;
    }

} // mpc