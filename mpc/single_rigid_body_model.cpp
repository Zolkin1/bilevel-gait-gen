//
// Created by zolkin on 1/18/24.
//

#include "pinocchio/algorithm/center-of-mass.hpp"
#include "pinocchio/algorithm/centroidal.hpp"

#include "single_rigid_body_model.h"
#include "timer.h"

namespace mpc {
    // TODO: For the SRBD we know the number of states at compile time so we should be able to remove all of the DMA

    SingleRigidBodyModel::SingleRigidBodyModel(const std::string& robot_urdf, const std::vector<std::string>& frames,
                                               int discretization_steps, double dt, const vector_t& nom_state) :
            Model(robot_urdf, frames, discretization_steps, dt, false,
                  {Constraints::Dynamics}),//, Constraints::ForceBox, Constraints::FrictionCone}), //, Constraints::EndEffectorLocation}),
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
            A.block<3,1>(ANG_VEL_START, ANG_VEL_START + i).noalias() =
                    -Id.col(i).cross(Ir_*omega) - omega.cross(Ir_.col(i));
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
                if (traj.IsForceMutable(ee, coord, time)) {
                    // TODO: DMA
                    vector_t vars_lin = traj.GetSplineLin(Trajectory::SplineTypes::Force, ee, coord, time);

                    int vars_idx, vars_affecting;
                    std::tie(vars_idx, vars_affecting) = traj.GetForceSplineIndex(ee, time, coord);

                    // linear momentum
                    B.block(LIN_MOM_START + coord, vars_idx - vars_affecting, 1, vars_affecting) =
                            vars_lin.transpose();

                    // Angular momentum - force
                    for (int poly = 0; poly < vars_lin.size(); poly++) {
                        B.block(ANG_VEL_START, vars_idx - vars_affecting + poly, 3, 1).noalias() =
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
                        B.block(ANG_VEL_START, pos_spline_start + vars_idx - vars_affecting + poly, 3, 1).noalias() =
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
        xdot.segment<POS_VARS>(LIN_MOM_START).noalias() = robot_mass_*GRAVITY; //vector_3t::Zero()

        // Orientation
        xdot.segment<3>(ORIENTATION_START) = Ir_inv_ * omega;

        // Angular momentum
        xdot.segment<3>(ANG_VEL_START) = vector_t::Zero(3);


        const vector_3t com_pos = GetCOMPosition(state_man);
        for (int i = 0; i < num_ee_; i++) {
            // COM linear momentum
            const vector_3t force = traj.GetForce(i, time);
            xdot.segment<POS_VARS>(LIN_MOM_START) += force;

            // Angular momentum
            xdot.segment<3>(ANG_VEL_START) += -omega.cross(Ir_*omega) +
                    (traj.GetEndEffectorLocation(i, time) - com_pos).cross(force);
        }

        return xdot;
    }

} // mpc