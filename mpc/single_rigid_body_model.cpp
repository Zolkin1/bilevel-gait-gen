//
// Created by zolkin on 1/18/24.
//

#include "pinocchio/algorithm/center-of-mass.hpp"

#include "single_rigid_body_model.h"
#include "timer.h"

namespace mpc {
    // TODO: For the SRBD we know the number of states at compile time so we should be able to remove all of the DMA

    SingleRigidBodyModel::SingleRigidBodyModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                                               int discretization_steps, double dt) :
            Model(robot_urdf, frames, discretization_steps, dt, false,
                  {Constraints::Dynamics, Constraints::ForceBox, Constraints::FrictionCone, Constraints::EndEffectorLocation}),
            num_tangent_states_(MOMENTUM_OFFSET + FLOATING_VEL_OFFSET),
            num_manifold_states_(MOMENTUM_OFFSET + FLOATING_BASE_OFFSET) {
        // Populate Ir
        // Calc Ir inv
    }

    void SingleRigidBodyModel::ConvertMPCStateToPinocchioState(const vector_t &state, Eigen::Ref<vector_t> q_pin) const {
        assert(state.size() == num_manifold_states_); // TODO: Check
        q_pin.head(POS_VARS) = state.segment<POS_VARS>(MOMENTUM_OFFSET);
        q_pin.segment<4>(POS_VARS) = state.segment<4>(MOMENTUM_OFFSET + POS_VARS);
    }

    Eigen::Vector3d SingleRigidBodyModel::GetCOMPosition(const mpc::vector_t& state) const {
        assert(state.size() == num_manifold_states_);
        return state.segment<3>(MOMENTUM_OFFSET);
    }

    void SingleRigidBodyModel::GetLinearDiscreteDynamics(const mpc::vector_t& state, const mpc::vector_t& ref_state,
                                                         const mpc::Trajectory& traj, double time, mpc::matrix_t& A,
                                                         mpc::matrix_t& B, mpc::vector_t& C) {
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



        // ------------------------------------------- //
        // --------------- Calculate A --------------- //
        // ------------------------------------------- //
        A_timer.StartTimer();
        Ac_.resize(num_tangent_states_, num_tangent_states_);
        Ac_.setZero();

        // COM position
        Ac_.block<3,3>(0, 3) = matrix_33t::Identity();

        // Orientation
        Ac_.block<3,3>(MOMENTUM_OFFSET, MOMENTUM_OFFSET + 3) = matrix_33t::Identity();

        matrix_33t Id = matrix_33t::Identity();
        for (int i = 0; i < 3; i++) {
            Ac_.block<3,1>(MOMENTUM_OFFSET + 3, MOMENTUM_OFFSET + 3 + i) =
                    Ir_inv_*(Id.col(i).cross(Ir_*omega) + omega.cross(Ir_.col(i)));
        }

        std::cout << Ac_ << std::endl;
//
//
//        Ac_.block(MOMENTUM_OFFSET, 0, POS_VARS, POS_VARS) = matrix_t::Identity(POS_VARS, POS_VARS)/robot_mass_;
//        Ac_.block(MOMENTUM_OFFSET + 3, 3, POS_VARS, POS_VARS) = matrix_t::Identity(POS_VARS, POS_VARS)/robot_mass_;
//        A_timer.StopTimer();
//
//        // ------------------------------------------- //
//        // --------------- Calculate B --------------- //
//        // ------------------------------------------- //
//        B_timer.StartTimer();
//        int num_inputs = input.GetNumInputs();
//
//        // Inputs: [f_{c1,p1}^{x}(1), ... , f_{c1, p1}^{x}(4), f_{c1,p2}^{x}(1), ... , f_{c1, p2}^{x}(4), ... ]
//        //          = [f_{c1,p1}^{x}, f_{c1,p2}^{x}, f_{c1,p3}^{x}, ... ]
//        //          = [f_{c1}^{x}, f_{c1}^{y}, f_{c1}^{z}, f_{c2}^{x}, f_{c2}^{y}, f_{c2}^{z}, ... ]
//        //     Given a specific time, we only need to look at one pj (polynomial j) and thus
//        //     one set of columns for each spline.
//
//        Bc_ = matrix_t::Zero(num_tangent_states_, num_inputs);
//        const Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();
//        for (int ee = 0; ee < num_ee_; ee++) {
//            const Eigen::Vector3d ee_pos_wrt_com = Eigen::Vector3d::Zero();//TODO: Get from spline
//            for (int coord = 0; coord < 3; coord++) {
//                if (input.IsForceMutable(ee, coord, time)) {
//                    vector_t vars_lin = input.GetForces().at(ee).at(coord).GetPolyVarsLin(time);
//
//                    int vars_idx, vars_affecting;
//                    std::tie(vars_idx, vars_affecting) = input.GetForceSplineIndex(ee, time, coord);
//
//                    // linear momentum
//                    Bc_.block(coord, vars_idx - vars_affecting, 1, vars_affecting) = vars_lin.transpose();
//
//                    // Angular momentum
//                    for (int poly = 0; poly < vars_lin.size(); poly++) {
//                        Bc_.block(3, vars_idx - vars_affecting + poly, 3, 1).noalias() =
//                                ee_pos_wrt_com.cross(static_cast<Eigen::Vector3d>(Id.col(coord))) * vars_lin(poly);
//                    }
//                }
//            }
//        }
//        B_timer.StopTimer();
//
//        // ------------------------------------------- //
//        // ---------------- Calculate C -------------- //
//        // ------------------------------------------- //
//        C_timer.StartTimer();
//        const vector_t state_alg = ConvertManifoldStateToAlgebraState(state, ref_state);
//        Cc_.noalias() = -Ac_*state_alg;
//        Cc2_ = Cc_;
//        Cc_.noalias() += CalcDynamics(state_alg, input, time, ref_state) - Bc_*input.AsQPVector(time);
//        Cc2_.noalias() += CalcDynamics(state_alg, input, time + integrator_->GetDt()/2, ref_state) - Bc_*input.AsQPVector(time+ integrator_->GetDt()/2);
//        C_timer.StopTimer();
//
//        // Discretize with the integrator
//        integrator_timer.StartTimer();
//        integrator_->CalcDerivWrtStateSingleStep(state, Ac_, A);
//        integrator_->CalcDerivWrtInputSingleStep(state, Bc_, Ac_, B);    // TODO: Technically this Ac should be evaluated at a different time
//        integrator_->CalcLinearTermDiscretization(Cc_, Cc2_, Ac_, C);
//        integrator_timer.StopTimer();
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

    vector_t SingleRigidBodyModel::CalcDynamics(const mpc::vector_t& state, const mpc::Trajectory& traj, double time,
                                                const mpc::vector_t& ref_state) {
        assert(state.size() == num_tangent_states_);
        assert(state.size() == ref_state.size());

        const vector_t state_man = ConvertTangentStateToManifoldState(state, ref_state);

        // COM position
        xdot_.head<3>() = state.segment<3>(3)/robot_mass_;

        // COM linear momentum
        xdot_.segment<3>(3).noalias() = robot_mass_*GRAVITY;

        // Orientation
        xdot_.segment<3>(6) = state.segment<3>(9);

        // Angular momentum
        xdot_.segment<3>(9) = vector_t::Zero(3);


        const Eigen::Vector3d com_pos = GetCOMPosition(state_man);
        for (int i = 0; i < num_ee_; i++) {
            // COM linear momentum
            const Eigen::Vector3d force = traj.GetForce(i, time);
            xdot_.topRows<3>() += force;

            // Angular momentum
            xdot_.segment<3>(6) += (traj.GetEndEffectorLocation(i, time) - com_pos).cross(force);
        }

        return xdot_;
    }

} // mpc