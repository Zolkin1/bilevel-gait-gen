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
#include "pinocchio/algorithm/center-of-mass.hpp"

#include "centroidal_model.h"
#include "euler_integrator.h"
#include "rk_integrator.h"

namespace mpc {
    using matrix6x_t = Eigen::Matrix<double, 6, Eigen::Dynamic>;

    CentroidalModel::CentroidalModel(const std::string &robot_urdf, const std::vector<std::string>& frames,
                                     int discretization_steps, double dt) :
            GRAVITY(0., 0., -9.81) {

        integrator_ = std::make_unique<RKIntegrator>(dt);

        // create the pinocchio model - always a free flyer
        pinocchio::urdf::buildModel(robot_urdf, pinocchio::JointModelFreeFlyer(), pin_model_, false);

        // create the pinocchio data
        pin_data_ = std::make_unique<pinocchio::Data>(pin_model_);

        num_joints_ = pin_model_.nq - FLOATING_BASE_OFFSET;
        num_total_states_ = MOMENTUM_OFFSET + FLOATING_VEL_OFFSET + num_joints_;

        robot_mass_ = pinocchio::computeTotalMass(pin_model_);

        num_ee_ = frames.size();
        frames_= frames;

        CreateFrameMap(frames);

        if (discretization_steps != 1) {
            throw std::runtime_error("Only discretization step of 1 is currently supported.");
        }
        discretization_steps_ = discretization_steps;

//        ref_state_ = vector_t::Zero(num_total_states_+1);
    }

    void CentroidalModel::GetLinearDiscreteDynamics(const vector_t& state, const vector_t& ref_state, const Inputs& input,
                                                             double time, matrix_t& A, matrix_t& B,
                                                             vector_t& C) {

        // Get the pinocchio configuration
        vector_t q_pin = ConvertMPCStateToPinocchioState(state);

        // ------------------------------------------- //
        // --------------- Calculate A --------------- //
        // ------------------------------------------- //
        matrix_t Ac = matrix_t::Zero(num_total_states_, num_total_states_);

        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q_pin);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q_pin);

        matrix_3t J = matrix_t::Zero(3, pin_model_.nv);
        for (int i = 0; i < num_ee_; i++) {
            // Get the jacobian of the FK
            J = GetFKJacobianForEndEffector(q_pin, frames_.at(i), false);
            Eigen::Vector3d force = input.GetForce(i, time);
            for (int j = 0; j < pin_model_.nv; j++) {
                // Cross product with the corresponding input force, add these for all the end effectors
                Ac.block<3,1>(POS_VARS, MOMENTUM_OFFSET + j) +=
                        static_cast<Eigen::Vector3d>(J.col(j)).cross(force);
            }
        }

        // Deal with centroidal dynamics
//        matrix6x_t CMM = pinocchio::computeCentroidalMap(pin_model_, *pin_data_, q_pin);
//
//        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
//        matrix_t CMMj = CMM.rightCols(num_joints_);
//        matrix_t CMMbinv = CMMb.inverse();  // TODO: maybe calc better, see below

        // OCS2 Version
//        const double mass = CMMb(0, 0);
//        Eigen::Matrix<double, 3, 3> Ab_22_inv = CMMb.template block<3, 3>(3, 3).inverse();
//        Eigen::Matrix<double, 6, 6> Ab_inv = Eigen::Matrix<double, 6, 6>::Zero();
//        Ab_inv << 1.0 / mass * Eigen::Matrix<double, 3, 3>::Identity(), -1.0 / mass * CMMb.template block<3, 3>(0, 3) * Ab_22_inv,
//                Eigen::Matrix<double, 3, 3>::Zero(), Ab_22_inv;
//
//        vector_t v_pin = ComputePinocchioVelocities(state, input.GetVels(time), CMMbinv, CMMj);

//        std::cout << "pinocchio velocities: \n" << v_pin << std::endl;

//        vector_t a_pin = vector_t::Zero(pin_model_.nv);
//
//        matrix6x_t dhdq_m = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
//        matrix6x_t dhdotdq = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
//        matrix6x_t dhdotdv = matrix6x_t::Zero(MOMENTUM_OFFSET, pin_model_.nv);
//
//        pinocchio::updateFramePlacements(pin_model_, *pin_data_);
        // Get dhdq
        // TODO: I only need this for dhdq so I can just steal that one line for efficiency, see below
//        pinocchio::computeCentroidalDynamicsDerivatives(pin_model_, *pin_data_, q_pin, v_pin, a_pin,
//                                                        dhdq_m, dhdotdq, dhdotdv, CMM);

        // OCS2 version
//        matrix6x_t dhdq(6, pin_model_.nv);
//        pinocchio::translateForceSet(pin_data_->dHdq, pin_data_->com[0], dhdq.const_cast_derived());
//        for (size_t k = 0; k < pin_model_.nv; ++k) {
//            dhdq.template block<3, 1>(pinocchio::Force::ANGULAR, k) +=
//                    pin_data_->hg.linear().cross(pin_data_->dFda.template block<3, 1>(pinocchio::Force::LINEAR, k))
//                    / pin_data_->Ig.mass();
//        }

//        std::cout << "dhdq me: \n" << dhdq_m << std::endl;
//        std::cout << "dhdq ocs2: \n" << dhdq << std::endl;

//        Ac.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) << CMMbinv, CMMbinv * dhdq;
        // TODO: Note CMM has been removed
        Ac.block(MOMENTUM_OFFSET, 0, POS_VARS, POS_VARS) = matrix_t::Identity(POS_VARS, POS_VARS)/robot_mass_;
        Ac.block(MOMENTUM_OFFSET + 3, 3, POS_VARS, POS_VARS) = matrix_t::Identity(POS_VARS, POS_VARS)/robot_mass_;

//        assert(std::abs(A(8, 2) - 1/robot_mass_) <= 1e-2);

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

        matrix_t Bc = matrix_t::Zero(num_total_states_, num_inputs);
        // --- Linear momentum --- //
        // Linear in each force, and each force is linear in its coefficients at each point in time.
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int coord = 0; coord < 3; coord++) {
                if (input.IsForceMutable(ee, coord, time)) {
                    vector_t vars_lin = input.GetForces().at(ee).at(coord).GetPolyVarsLin(time);

                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) = input.GetForceSplineIndex(ee, time, coord);

                    Bc.block(coord, vars_idx - vars_affecting, 1, vars_affecting) = vars_lin.transpose();
                }
            }
        }

        const Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();
        // --- Angular momentum --- //
        vector_t vars_lin;
        for (int ee = 0; ee < num_ee_; ee++) {
            Eigen::Vector3d ee_pos_wrt_com = GetEndEffectorLocationCOMFrame(state, frames_.at(ee));
            for (int coord = 0; coord < 3; coord++) {
                if (input.IsForceMutable(ee, coord, time)) {
                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) = input.GetForceSplineIndex(ee, time, coord);
                    vars_lin = input.GetForces().at(ee).at(coord).GetPolyVarsLin(time);

                    for (int poly = 0; poly < vars_lin.size(); poly++) {
                        Bc.block(3, vars_idx - vars_affecting + poly, 3, 1) =
                                ee_pos_wrt_com.cross(static_cast<Eigen::Vector3d>(Id.col(coord))) * vars_lin(poly);
                    }
                }

            }
        }

//        std::cout << "B: \n" << B << std::endl;

        // --- Base velocity --- //
//        Bc.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) <<
//                matrix_t::Zero(FLOATING_VEL_OFFSET, num_inputs - num_joints_),
//                -CMMbinv*CMMj;

        // --- Joint velocity --- //
        Bc.bottomRows(num_joints_) << matrix_t::Zero(num_joints_, num_inputs - num_joints_),
                matrix_t::Identity(num_joints_, num_joints_);


        // ------------------------------------------- //
        // ---------------- Calculate C -------------- //
        // ------------------------------------------- //
        vector_t state_alg = ConvertManifoldStateToAlgebraState(state, ref_state);
        // TODO: move this to an integrator function
        vector_t Cc = CalcDynamics(state_alg, input, time, ref_state) -Ac*state_alg -Bc*input.AsQPVector(time);
        vector_t Cc2 = CalcDynamics(state_alg, input, time + integrator_->GetDt()/2, ref_state)
                -Ac*state_alg -Bc*input.AsQPVector(time+ integrator_->GetDt()/2);

        // Discretize with the integrator
        A = integrator_->CalcDerivWrtStateSingleStep(state, Ac);
        B = integrator_->CalcDerivWrtInputSingleStep(state, Bc, Ac);    // TODO: Technically this Ac should be evaluated at a different time
        C = integrator_->CalcLinearTermDiscretization(Cc, Cc2, Ac); //GetDt()*Cc;

        // TODO: Remove
        for (int i = 0; i < C.size(); i++) {
            assert(!std::isnan(C(i)));
        }
    }

    void CentroidalModel::GetFKLinearization(const vector_t& state, const vector_t& ref_state, const Inputs& input, int end_effector,
                                                 matrix_t& A, vector_t& C) {
        const vector_t q_pin = ConvertMPCStateToPinocchioState(state);
        A = matrix_t::Zero(3, pin_model_.nv);
        C = vector_t::Zero(3);

        pinocchio::computeJointJacobians(pin_model_, *pin_data_, q_pin);
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q_pin);

        matrix6x_t J = matrix6x_t::Zero(6, pin_model_.nv);
        pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_map_.at(frames_.at(end_effector)),
                                    pinocchio::LOCAL_WORLD_ALIGNED, J);

        A = J.topRows<3>();// - Eigen::Matrix3Xd::Identity(3, pin_model_.nv);

        const vector_t state_alg = ConvertManifoldStateToAlgebraState(state, ref_state);
        C = GetEndEffectorLocationCOMFrame(state, frames_.at(end_effector)) + GetCOMPosition(state)
                - A*state_alg.tail(pin_model_.nv);
    }

    matrix_3t CentroidalModel::GetFKJacobianForEndEffector(const vector_t& q, const std::string& frame, bool compute_jac) {
        if (compute_jac) {
            pinocchio::computeJointJacobians(pin_model_, *pin_data_, q);
            pinocchio::framesForwardKinematics(pin_model_, *pin_data_, q);
        }

        matrix6x_t J = matrix_t::Zero(6,pin_model_.nv);
        pinocchio::getFrameJacobian(pin_model_, *pin_data_, frame_map_.at(frame), pinocchio::LOCAL_WORLD_ALIGNED, J);

        // Only return the position coordinates
        return J.topRows<3>() - Eigen::Matrix3Xd::Identity(3, pin_model_.nv);
    }

    Eigen::Vector3d CentroidalModel::GetEndEffectorLocationCOMFrame(const vector_t& state, const std::string& frame) const {
        pinocchio::forwardKinematics(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state));
        pinocchio::framesForwardKinematics(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state));
        return pin_data_->oMf.at(frame_map_.at(frame)).translation() - GetCOMPosition(state);
    }

    vector_t CentroidalModel::ComputeBaseVelocities(const vector_t& state, const vector_t& vel) const {
        const vector_t q_pin = ConvertMPCStateToPinocchioState(state);
        const matrix6x_t CMM = pinocchio::computeCentroidalMap(pin_model_, *pin_data_, q_pin);

        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
        matrix_t CMMj = CMM.rightCols(num_joints_);
        matrix_t CMMbinv = CMMb.inverse();  // TODO: maybe calc better
        return ComputePinocchioVelocities(state, vel, CMMbinv, CMMj);
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

        v_pin.head<FLOATING_VEL_OFFSET>().noalias() = CMMbinv * (state.head<MOMENTUM_OFFSET>() - CMMj*joint_vels);

        v_pin.tail(num_joints_) = joint_vels;

        return v_pin;
    }

    vector_t CentroidalModel::ConvertMPCStateToPinocchioState(const vector_t &state) const {
        assert(state.size() == pin_model_.nq + 6);
        vector_t q(pin_model_.nq);

        q.head(POS_VARS) = state.segment<POS_VARS>(MOMENTUM_OFFSET);
        q.segment<4>(POS_VARS) = state.segment<4>(MOMENTUM_OFFSET + POS_VARS);
        q.tail(num_joints_) = state.tail(num_joints_);

        return q;
    }

    vector_t CentroidalModel::ConvertPinocchioStateToMPCState(const Eigen::Vector<double, MOMENTUM_OFFSET>& momentum,
                                                              const vector_t& state) const {
        const Eigen::Vector3d rot = ConvertQuaternionToZYXRot(state.segment<4>(POS_VARS));
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
        assert(state.size() == ref_state.size());
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
        assert(state.size() == ref_state.size()-1);
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

    Eigen::Vector3d CentroidalModel::GetCOMPosition(const vector_t& state) const {
        assert(state.size() == pin_model_.nq + MOMENTUM_OFFSET);
//        pinocchio::forwardKinematics(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state));
//        return pin_data_->com[0];
        return pinocchio::centerOfMass(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state));
    }

    void CentroidalModel::CreateFrameMap(const std::vector<std::string>& frames) {
        for (int i = 0; i < num_ee_; i++) {
            for (int j = 0; j < pin_model_.frames.size(); j++)
            if (pin_model_.frames.at(j).name == frames.at(i)) {
                frame_map_.insert(std::pair<std::string, int>(frames.at(i), pin_model_.getFrameId(frames.at(i))));
            }
        }
    }

    vector_t CentroidalModel::CalcDynamics(const vector_t& state, const Inputs& input, double time,
                                           const vector_t& ref_state) const {
        assert(state.size() == pin_model_.nv + 6);

        const vector_t state_man = ConvertAlgebraStateToManifoldState(state, ref_state);

        vector_t xdot = vector_t::Zero(num_total_states_);

        xdot.topRows<3>() = robot_mass_*GRAVITY;

        for (int i = 0; i < num_ee_; i++) {
            xdot.topRows<3>() += input.GetForce(i, time);

            // Note that the input ee position can deviate from the state ee position.
            // So I use the state ee position.
            xdot.middleRows<3>(3) += GetEndEffectorLocationCOMFrame(state_man, frames_.at(i)).cross(input.GetForce(i, time));
        }

        // NOTE: This CMM call is expensive
//        matrix_t CMM = pinocchio::computeCentroidalMap(pin_model_, *pin_data_, ConvertMPCStateToPinocchioState(state_man));
//        matrix_t CMMb = CMM.leftCols<FLOATING_VEL_OFFSET>();
//        matrix_t CMMj = CMM.rightCols(num_joints_);
//        xdot.middleRows<FLOATING_VEL_OFFSET>(MOMENTUM_OFFSET) =
//                CMMb.householderQr().solve(state.head<MOMENTUM_OFFSET>() - CMMj*input.GetVels(time));
        xdot.segment<POS_VARS>(MOMENTUM_OFFSET) = state.segment<POS_VARS>(0)/robot_mass_;
        xdot.segment<POS_VARS>(MOMENTUM_OFFSET + POS_VARS) = state.segment<POS_VARS>(POS_VARS)/robot_mass_;

        xdot.bottomRows(num_joints_) = input.GetVels(time);

        return xdot;
    }

    double CentroidalModel::GetMass() const {
        return robot_mass_;
    }

    vector_t CentroidalModel::GetDiscreteDynamics(const vector_t& state, const Inputs& input, double time,
                                                  const vector_t& ref_state) const {
        vector_t xkp1 = integrator_->CalcIntegral(state, input, time, 1, *this, ref_state); //state + integrator_->GetDt()*CalcDynamics(state, input, time, ref_state);
        return xkp1;
    }

    const std::string& CentroidalModel::GetEndEffectorFrame(int ee) const {
        return frames_.at(ee);
    }

    std::vector<int> CentroidalModel::GetContactFrames() const {
        std::vector<int> frames;
        for (const auto& frame : frames_) {
            frames.push_back(frame_map_.at(frame));
        }

        return frames;
    }

} // mpc