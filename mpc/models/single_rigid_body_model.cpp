//
// Created by zolkin on 1/18/24.
//

#include "pinocchio/algorithm/center-of-mass.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"

#include "models/single_rigid_body_model.h"
#include "timer.h"

namespace mpc {
    // TODO: For the SRBD we know the number of states at compile time so we should be able to remove all of the DMA

    SingleRigidBodyModel::SingleRigidBodyModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                                               int discretization_steps, double dt, const vector_t& nom_state) :
            Model(robot_urdf, frames, discretization_steps, dt, false,
                  {Constraints::Dynamics,
                   Constraints::ForceBox,
                   Constraints::FrictionCone,
                   Constraints::EndEffectorLocation
                  }),
            num_tangent_states_(MOMENTUM_OFFSET + FLOATING_VEL_OFFSET),
            num_manifold_states_(MOMENTUM_OFFSET + FLOATING_BASE_OFFSET) {
        // Populate Ir
        computeCentroidalMap(pin_model_, *pin_data_, nom_state);
        Ir_ = pin_data_->oMi[1].actInv(pin_data_->oYcrb[0]).inertia();

        // Calc Ir inv
        Ir_inv_ = Ir_.inverse();
    }

    void SingleRigidBodyModel::ConvertMPCStateToPinocchioState(const vector_t &state, Eigen::Ref<vector_t> q_pin) const {
        assert(state.size() == num_manifold_states_); // TODO: Check
        q_pin.head(POS_VARS) = state.segment<POS_VARS>(POS_START);
        q_pin.segment<4>(POS_VARS) = state.segment<4>(QUAT_START);
    }

    Eigen::Vector3d SingleRigidBodyModel::GetCOMPosition(const mpc::vector_t& state) const {
        assert(state.size() == num_manifold_states_);
        return state.segment<3>(POS_START);
    }

    void SingleRigidBodyModel::GetLinearDynamics(const vector_t& state,
                                                 const vector_t& ref_state,
                                                 const Trajectory& traj,
                                                 double dt,
                                                 double time,
                                                 matrix_t& A, matrix_t& B, vector_t& C, vector_t& C2) {
        utils::Timer A_timer("A creation");
        utils::Timer B_timer("B creation");
        utils::Timer C_timer("C creation");
        utils::Timer integrator_timer("integration");

        assert(state.size() == num_manifold_states_);
        // x: [p, pdot, dq, omega]
        const vector_3t& p = state.head<3>();
        const vector_3t& pdot = state.segment<3>(3);
        const vector_3t dq = ConvertManifoldToTangentQuat(state.segment<QUAT_SIZE>(QUAT_START), ref_state.segment<QUAT_SIZE>(QUAT_START));
        const vector_3t& omega = state.tail<3>();

        const matrix_33t Id = matrix_33t::Identity();

        // ------------------------------------------- //
        // --------------- Calculate A --------------- //
        // ------------------------------------------- //
        A_timer.StartTimer();
        A.resize(num_tangent_states_, num_tangent_states_);
        A.setZero();

        // COM position
        A.block<3,3>(POS_START, LIN_MOM_START) = matrix_33t::Identity()/robot_mass_;

        // Orientation
        A.block<3,3>(ORIENTATION_START, ANG_VEL_START) = Ir_inv_;

        // Angular momentum
        for (int i = 0; i < 3; i++) {
            // Omega dependence
            A.block<3,1>(ANG_VEL_START, ANG_VEL_START + i).noalias() =
                    -Id.col(i).cross(Ir_*omega) - omega.cross(Ir_.col(i));

            // COM dependence
            for (int ee = 0; ee < num_ee_; ee++) {
                A.block<3, 1>(ANG_VEL_START, POS_START + i) += -Id.col(i).cross(traj.GetForce(ee, time));
            }
        }

//        std::cout << A << std::endl;
        A_timer.StopTimer();

        // ------------------------------------------- //
        // --------------- Calculate B --------------- //
        // ------------------------------------------- //
        B_timer.StartTimer();
        const int num_inputs = traj.GetTotalForceSplineVars() + traj.GetTotalPosSplineVars();
        const int pos_spline_start = traj.GetTotalForceSplineVars();

        B.resize(num_tangent_states_, num_inputs);
        B.setZero();

        for (int ee = 0; ee < num_ee_; ee++) {
            const vector_3t ee_pos_wrt_com = traj.GetEndEffectorLocation(ee, time) - GetCOMPosition(state);
            const vector_3t force = traj.GetForce(ee, time);

            for (int coord = 0; coord < 3; coord++) {
                if (traj.IsForceMutable(ee, time)) {
                    // TODO: DMA
                    vector_t vars_lin = traj.GetSplineLin(Trajectory::SplineTypes::Force, ee, coord, time);

                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) = traj.GetForceSplineIndex(ee, time, coord);

                    // linear momentum
                    B.block(LIN_MOM_START + coord, vars_idx, 1, vars_affecting) =
                            vars_lin.transpose();

                    // Angular momentum - force
                    for (int poly = 0; poly < vars_lin.size(); poly++) {
                        B.block(ANG_VEL_START, vars_idx + poly, 3, 1).noalias() =
                                ee_pos_wrt_com.cross(static_cast<Eigen::Vector3d>(Id.col(coord))) * vars_lin(poly);
                    }
                }

                // Angular momentum - position
                if (coord != 2) {
                    // TODO: DMA
                    vector_t vars_lin = traj.GetSplineLin(Trajectory::SplineTypes::Position, ee, coord, time);
                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) = traj.GetPositionSplineIndex(ee, time, coord);
                    for (int poly = 0; poly < vars_lin.size(); poly++) {
                        B.block(ANG_VEL_START, pos_spline_start + vars_idx + poly, 3, 1).noalias() =
                                Id.col(coord).cross(force) * vars_lin(poly);
                    }
                }
            }
        }
        B_timer.StopTimer();
//        std::cout << B << std::endl;


        // ------------------------------------------- //
        // ---------------- Calculate C -------------- //
        // ------------------------------------------- //
        C_timer.StartTimer();
        const vector_t tan_state = ConvertManifoldStateToTangentState(state, ref_state);
        C.noalias() = -A*tan_state - B*traj.SplinesAsVec();
        C2 = C;
        C += CalcDynamics(tan_state, traj, time, ref_state);
        C2 += CalcDynamics(tan_state, traj, time + dt/2, ref_state);
        C_timer.StopTimer();

//        vector_t xdot = CalcDynamics(tan_state, traj, time, ref_state);
//        vector_t lin_dyn = A*tan_state + B*traj.SplinesAsVec() + C;
//
//        std::cout << "diff: \n" << xdot - lin_dyn << std::endl;

    }

    int SingleRigidBodyModel::GetNumTangentStates() const {
        return num_tangent_states_;
    }

    int SingleRigidBodyModel::GetNumManifoldStates() const {
        return num_manifold_states_;
    }

    vector_3t SingleRigidBodyModel::ConvertManifoldToTangentQuat(const Eigen::Vector4d& state,
                                                                 const Eigen::Vector4d& ref_state) {
        Eigen::Quaterniond quat(state);
        Eigen::Quaterniond quat_ref(ref_state);
        vector_3t tangent_state = pinocchio::quaternion::log3(quat_ref.inverse()*quat);

        return tangent_state;
    }

    vector_t SingleRigidBodyModel::ConvertManifoldStateToTangentState(const mpc::vector_t& state,
                                                                      const mpc::vector_t& ref_state) const {
        assert(state.size() == ref_state.size());
        // TODO: DMA
        vector_t tangent_state(state.size()-1);
        tangent_state.head<QUAT_START>() = state.head<QUAT_START>();
        tangent_state.segment<QUAT_SIZE-1>(QUAT_START) =
                ConvertManifoldToTangentQuat(state.segment<QUAT_SIZE>(QUAT_START),
                        ref_state.segment<QUAT_SIZE>(QUAT_START));
        tangent_state.tail<3>() = state.tail<3>();

        return tangent_state;
    }

    vector_t SingleRigidBodyModel::ConvertTangentStateToManifoldState(const mpc::vector_t& state,
                                                                      const mpc::vector_t& ref_state) const {
        // TODO: DMA
        vector_t manifold_state(state.size() + 1);
        manifold_state.head<QUAT_START>() = state.head<QUAT_START>();

        Eigen::Quaterniond quat;
        pinocchio::quaternion::exp3(state.segment<3>(QUAT_START), quat);
        Eigen::Quaterniond quat_ref(static_cast<Eigen::Vector4d>(ref_state.segment<4>(QUAT_START)));
        quat = quat_ref*quat;
        manifold_state(QUAT_START) = quat.x();
        manifold_state(1+QUAT_START) = quat.y();
        manifold_state(2+QUAT_START) = quat.z();
        manifold_state(3+QUAT_START) = quat.w();

        manifold_state.tail<3>() = state.tail<3>();

        return manifold_state;
    }

    tan_state_t SingleRigidBodyModel::CalcDynamics(const mpc::vector_t& state, const mpc::Trajectory& traj, double time,
                                                   const mpc::vector_t& ref_state) {
        assert(state.size() == num_tangent_states_);
        assert(state.size() == ref_state.size() - 1);

        const vector_t state_man = ConvertTangentStateToManifoldState(state, ref_state);
        const vector_3t& omega = state.segment<3>(ANG_VEL_START);

        tan_state_t xdot;

        // COM position
        xdot.head<POS_VARS>() = state.segment<POS_VARS>(LIN_MOM_START)/robot_mass_;

        // COM linear momentum
        xdot.segment<POS_VARS>(LIN_MOM_START).noalias() = robot_mass_*GRAVITY; //vector_3t::Zero();

        // Orientation
        xdot.segment<3>(ORIENTATION_START).noalias() = Ir_inv_ * omega;

        // Angular momentum
        xdot.segment<3>(ANG_VEL_START).noalias() = -omega.cross(Ir_*omega);


        const vector_3t com_pos = GetCOMPosition(state_man);
        for (int i = 0; i < num_ee_; i++) {
            // COM linear momentum
            const vector_3t force = traj.GetForce(i, time);
            xdot.segment<POS_VARS>(LIN_MOM_START) += force;

            // Angular momentum
            xdot.segment<3>(ANG_VEL_START) += (traj.GetEndEffectorLocation(i, time) - com_pos).cross(force);
        }

        return xdot;
    }

    vector_3t SingleRigidBodyModel::GetCOMToHip(int end_effector) const {
        // TODO: Make not hard coded

//        for (int i = 0; i < pin_model_.njoints; i++) {
//            std::cout << pin_model_.names.at(i) << std::endl;
//        }

        int com_joint_id = pin_model_.getJointId("root_joint");

        // TODO: Consider getting a different frame/joint
        int hip_joint_id;
        switch (end_effector) {
            case 0:
                hip_joint_id = pin_model_.getJointId("FL_hip_joint");
                break;
            case 1:
                hip_joint_id = pin_model_.getJointId("FR_hip_joint");
                break;
            case 2:
                hip_joint_id = pin_model_.getJointId("RL_hip_joint");
                break;
            case 3:
                hip_joint_id = pin_model_.getJointId("RR_hip_joint");
                break;
            default:
                throw std::runtime_error("Not a valid end effector number");
        }

//        std::cout << "hip (wrt com): \n " << pin_data_->oMi[com_joint_id].translation()
//                                            - pin_data_->oMi[hip_joint_id].translation();

// TODO: This is adding an additional shift
        vector_3t temp = -pin_data_->oMi[com_joint_id].translation() + pin_data_->oMi[hip_joint_id].translation();
        if (temp(1) >= 0) {
            temp(1) += 0.13;
        } else {
            temp(1) -= 0.13;
        }

        return temp;
    }

    vector_t SingleRigidBodyModel::InverseKinematics(const mpc::man_state_t& state,
                                                     const std::vector<vector_3t>& end_effector_location,
                                                     const vector_t& state_guess,
                                                     const vector_t& joint_limits_ub,
                                                     const vector_t& joint_limits_lb) {
        utils::Timer ik_timer("IK solve");
        ik_timer.StartTimer();

        const int base_joint_id = pin_model_.getJointId("root_joint");
        const int num_constraints = 6+(3*num_ee_);

        assert(end_effector_location.size() == num_ee_);

        std::vector<int> frames;
        std::vector<pinocchio::SE3> ee_des;
        for (int ee = 0; ee < num_ee_; ee++) {
            frames.push_back(frame_map_.at(frames_.at(ee)));
            ee_des.emplace_back(matrix_33t::Identity(), end_effector_location.at(ee));
        }
//        const pinocchio::SE3 ee_des(Eigen::Matrix3d::Identity(), end_effector_location);
        const pinocchio::SE3 body_des(static_cast<Eigen::Quaterniond>(state.segment<QUAT_SIZE>(QUAT_START)).matrix(),
                                      state.segment<POS_VARS>(POS_START));

        Eigen::VectorXd q = pinocchio::neutral(pin_model_);
        q.head<POS_VARS>() = state.head<POS_VARS>();
        q.segment<QUAT_SIZE>(POS_VARS) = state.segment<QUAT_SIZE>(QUAT_START);
        q.tail(pin_model_.njoints-2) = state_guess.tail(pin_model_.njoints-2);

        const double eps  = 1e-4;
        const int IT_MAX  = 1000;
        const double DT   = 1e-1;
        const double damp = 1e-6;

        pinocchio::Data::Matrix6x Jee(6,pin_model_.nv);
        Jee.setZero();

        Eigen::Matrix<double, 6, Eigen::Dynamic> Jb(6, pin_model_.nv);
        Jb.setZero();

        matrix_t J(9, pin_model_.nv);
        J.setZero();

        bool success = false;
        typedef Eigen::Matrix<double, 6, 1> Vector6d;
        vector_t err(9);

        Eigen::VectorXd v(pin_model_.nv);

        std::vector<pinocchio::SE3> ee_err(num_ee_);

        // TODO: Figure out how to do this in one IK solve
        for (int ee = 0; ee < num_ee_; ee++) {
            for (int i = 0; i < IT_MAX; i++) {

                Eigen::Quaterniond quat(q.segment<4>(3));
                pinocchio::quaternion::firstOrderNormalize(quat);
                q(3) = quat.x();
                q(3 + 1) = quat.y();
                q(3 + 2) = quat.z();
                q(3 + 3) = quat.w();

                pinocchio::forwardKinematics(pin_model_, *pin_data_, q);
                pinocchio::updateFramePlacements(pin_model_, *pin_data_);

                ee_err.at(ee) = pin_data_->oMf[frames.at(ee)].actInv(
                        ee_des.at(ee)); //GetSE3Error(frame_id, ee_des);
                err.segment<3>(0) = ee_err.at(ee).translation();

                const pinocchio::SE3 body_err = pin_data_->oMi[base_joint_id].actInv(
                        body_des); //GetSE3Error(base_joint_id, body_des);
                err.tail<6>() = pinocchio::log6(body_err).toVector();

                if (err.norm() < eps) {
                    success = true;
                    break;
                }

                ComputeJacobianForIK(q, ee_err.at(ee), frames.at(ee), Jee, true);
                J.block(0, 0, 3, pin_model_.nv) = -Jee.topRows<3>();
                Jee.setZero();

                ComputeJacobianForIK(q, body_err, base_joint_id, Jb, false);
                J.block(3, 0, 6, pin_model_.nv) = Jb;
                Jb.setZero();

                matrix_t JJt(9, 9);
                JJt.noalias() = J * J.transpose();
                JJt.diagonal().array() += damp;
                v.noalias() = -(J.transpose() * JJt.ldlt().solve(err));
                q = pinocchio::integrate(pin_model_, q, v * DT);
                // TODO: Put in joint bounds
//                for (int j = 0; j < q.size() - 7; j++) {
//                    if (q(7+j) <= joint_limits_lb(j)) {
//                        q(7+j) = joint_limits_lb(j);
//                    } else if(q(7+j) >= joint_limits_ub(j)) {
//                        q(7+j) = joint_limits_ub(j);
//                    }
//                }
            }

            if (!success) {
                std::cerr << "IK did not converge." << std::endl;
            }
        }

        ik_timer.StopTimer();
//        ik_timer.PrintElapsedTime();

//        return state_out;
        return q;
    }

    pinocchio::SE3 SingleRigidBodyModel::GetSE3Error(int joint_id, const pinocchio::SE3& des_state) {
        return pin_data_->oMi[joint_id].actInv(des_state);
    }

    void SingleRigidBodyModel::ComputeJacobianForIK(const mpc::vector_t& q, const pinocchio::SE3& error,
                                                    int id, Eigen::Matrix<double, 6, Eigen::Dynamic>& J,
                                                    bool is_frame) {
        if (is_frame) {
            pinocchio::computeFrameJacobian(pin_model_, *pin_data_, q, id, J);
        } else {
            pinocchio::computeJointJacobian(pin_model_, *pin_data_, q, id, J);
            pinocchio::Data::Matrix6 Jlog;
            pinocchio::Jlog6(error.inverse(), Jlog);
            J = -Jlog * J;
        }
    }

    std::vector<vector_3t> SingleRigidBodyModel::GetEndEffectorLocations(const vector_t& q) {
        std::vector<vector_3t> ee_locations(4);
        assert(q.size() == pin_model_.nq);

        pinocchio::forwardKinematics(pin_model_, *pin_data_, q);
        pinocchio::updateFramePlacements(pin_model_, *pin_data_);

        for (int ee = 0; ee < num_ee_; ee++) {
            ee_locations.at(ee) = pin_data_->oMf[frame_map_.at(frames_.at(ee))].translation();
        }

        return ee_locations;
    }

    void SingleRigidBodyModel::ComputeLinearizationPartialWrtContactTimes(matrix_t& dA, matrix_t& dB,
                                                                 vector_t& dC, const vector_t& state,
                                                                 const Trajectory& traj, double time,
                                                                 int end_effector, int contact_time_idx) {

        // TODO: Can do this quicker by using the stuff computed for the original linearization

        const int num_inputs = traj.GetTotalPosSplineVars() + traj.GetTotalForceSplineVars();

        dA.resize(num_tangent_states_, num_tangent_states_);
        dB.resize(num_tangent_states_, num_inputs);
        dC.resize(num_tangent_states_);

        dA.setZero();
        dB.setZero();
        dC.setZero();

        const matrix_33t Id = matrix_33t::Identity();
        const vector_3t force_partial = traj.GetForcePartialWrtContactTime(end_effector, time, contact_time_idx);
        const vector_3t position_partial = traj.GetPositionPartialWrtContactTime(end_effector, time, contact_time_idx);

        // dA
        // COM position -- no dependence on contact times
        // Linear momentum -- no dependence in A
        // Orientation -- no dependence on contact times
        // Angular momentum
        for (int coord = 0; coord < POS_VARS; coord++) {
            // Omega dependence -- no dependence on contact times
            // TODO: confirm I shouldn't use the chain rule to back out to the contact times

            // COM dependence -- force appears here, which depends on contact times
            dA.block<3, 1>(ANG_VEL_START, POS_START + coord) += -Id.col(coord).cross(force_partial);
        }

        // dB
        const vector_3t ee_pos_wrt_com = traj.GetEndEffectorLocation(end_effector, time) - GetCOMPosition(state);
        const vector_3t force = traj.GetForce(end_effector, time);
        const int pos_spline_start = traj.GetTotalForceSplineVars();

        for (int coord = 0; coord < 3; coord++) {
            if (traj.IsForceMutable(end_effector, time)) {
                // TODO: DMA
                const vector_t force_coef_partials =
                        traj.GetForceCoefPartialsWrtContactTime(end_effector, coord, time, contact_time_idx);

                const vector_t vars_lin = traj.GetSplineLin(Trajectory::SplineTypes::Force, end_effector, coord, time);

                int vars_idx, vars_affecting;
                std::tie(vars_idx, vars_affecting) = traj.GetForceSplineIndex(end_effector, time, coord);

                // linear momentum -- force coefficients depend on contact times
                dB.block(LIN_MOM_START + coord, vars_idx, 1, vars_affecting) =
                        force_coef_partials.transpose(); // Not agreeing with the finite difference

                // Angular momentum - force -- force coefficients depend on contact times
                for (int poly = 0; poly < force_coef_partials.size(); poly++) {
                    dB.block(ANG_VEL_START, vars_idx + poly, 3, 1).noalias() =
                            ee_pos_wrt_com.cross(static_cast<Eigen::Vector3d>(Id.col(coord))) * force_coef_partials(poly) +
                            position_partial.cross(static_cast<Eigen::Vector3d>(Id.col(coord))) * vars_lin(poly);
                }
            }

            // Angular momentum - position -- position coefficients depend on contact times
            if (coord != 2) {
                // TODO: DMA
                const vector_t pos_coef_partials =
                        traj.GetPositionCoefPartialsWrtContactTime(end_effector, coord, time, contact_time_idx);
                const vector_t vars_lin = traj.GetSplineLin(Trajectory::SplineTypes::Position, end_effector, coord, time);

                int vars_idx, vars_affecting;
                std::tie(vars_idx, vars_affecting) = traj.GetPositionSplineIndex(end_effector, time, coord);
                for (int poly = 0; poly < pos_coef_partials.size(); poly++) {
                    dB.block(ANG_VEL_START, pos_spline_start + vars_idx + poly, 3, 1).noalias() =
                            Id.col(coord).cross(force) * pos_coef_partials(poly) +
                            Id.col(coord).cross(force_partial) * vars_lin(poly);
                }
            }
        }

        // TODO: Do this more efficiently
//        matrix_t A, B;
//        vector_t C, C2;
//        GetLinearDynamics(state, traj.GetState(0), traj, 0.02, time, A, B, C, C2);  // TODO: Make the dt not hard coded!

        // dC
        tan_state_t tan_state = ConvertManifoldStateToTangentState(state, traj.GetState(0));
        dC.noalias() = -dA*tan_state - dB*traj.SplinesAsVec();

        // Linear momentum depends on the force
        dC.segment<POS_VARS>(LIN_MOM_START) += force_partial;

        // Angular momentum depends on force and position
//        for (int coord = 0; coord < POS_VARS; coord++) {
            // position and force
        dC.segment<POS_VARS>(ANG_VEL_START).noalias() +=
                    ee_pos_wrt_com.cross(force_partial) + position_partial.cross(force);
//        }
    }

    SingleRigidBodyModel& SingleRigidBodyModel::operator=(const mpc::SingleRigidBodyModel& model) {
        if (this != &model) {
            num_tangent_states_ = model.num_tangent_states_;
            num_manifold_states_ = model.num_manifold_states_;
            Ir_ = model.Ir_;
            Ir_inv_ = model.Ir_inv_;
        }

        return *this;
    }

    SingleRigidBodyModel::SingleRigidBodyModel(const mpc::SingleRigidBodyModel& other) : Model(other) {
        *this = other;
    }


//    vector_t SingleRigidBodyModel::VelocityInverseKinematics(const vector_t& config,
//                                                             const std::vector<vector_3t>& end_effector_velocity,
//                                                             const vector_t& vel_guess) {
//
//        std::vector<int> frames;
//        for (int ee = 0; ee < num_ee_; ee++) {traj.GetForcePartialWrtContactTime(en
//            frames.push_back(frame_map_.at(frames_.at(ee)));
//            ee_des.emplace_back(matrix_33t::Identity(), end_effector_location.at(ee));
//        }
//
//        const double eps  = 1e-4;
//        const int IT_MAX  = 1000;
//        const double DT   = 1e-1;
//        const double damp = 1e-6;
//
//        pinocchio::Data::Matrix6x Jdotee(6,pin_model_.nv);
//        Jdotee.setZero();
//
//        Eigen::Matrix<double, 6, Eigen::Dynamic> Jdotb(6, pin_model_.nv);
//        Jdotb.setZero();
//
//        matrix_t Jdot(9, pin_model_.nv);
//        Jdot.setZero();
//
//        bool success = false;
//        typedef Eigen::Matrix<double, 6, 1> Vector6d;
//        vector_t err(9);
//
//        Eigen::VectorXd v(pin_model_.nv);
//        v.setZero();
//
//        Eigen::VectorXd a(pin_model_.nv);
//        a.setZero();
//
//        std::vector<vector_3t> ee_err(num_ee_);
//
//        // TODO: Figure out how to do this in one IK solve
//        for (int ee = 0; ee < num_ee_; ee++) {
//            for (int i = 0; i < IT_MAX; i++) {
//                pinocchio::forwardKinematics(pin_model_, *pin_data_, config, v);
//                pinocchio::updateFramePlacements(pin_model_, *pin_data_);
//
//                ee_err.at(ee) = pin_data_->oMf[frames.at(ee)].translation() - J*qdot;
//                err.segment<3>(0) = ee_err.at(ee);
//
//                const pinocchio::SE3 body_err = pin_data_->oMi[base_joint_id].actInv(
//                        body_des); //GetSE3Error(base_joint_id, body_des);
//                err.tail<6>() = pinocchio::log6(body_err).toVector();
//
//                if (err.norm() < eps) {
//                    success = true;
//                    break;
//                }
//
//                ComputeJacobianForIK(q, ee_err.at(ee), frames.at(ee), Jee, true);
//                J.block(0, 0, 3, pin_model_.nv) = -Jee.topRows<3>();
//                Jee.setZero();
//
//                ComputeJacobianForIK(q, body_err, base_joint_id, Jb, false);
//                J.block(3, 0, 6, pin_model_.nv) = Jb;
//                Jb.setZero();
//
//                matrix_t JJt(9, 9);
//                JJt.noalias() = J * J.transpose();
//                JJt.diagonal().array() += damp;
//                v.noalias() = -(J.transpose() * JJt.ldlt().solve(err));
//                q = pinocchio::integrate(pin_model_, q, v * DT);
//            }
//
//            if (!success) {
//                std::cerr << "IK did not converge." << std::endl;
//            }
//        }
//
//        ik_timer.StopTimer();
//
//        return q;
//    }

} // mpc