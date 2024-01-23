//
// Created by zolkin on 12/8/23.
//

#include "rk_integrator.h"
#include "centroidal_model.h"
#include "single_rigid_body_model.h"

namespace mpc {
    RKIntegrator::RKIntegrator(double dt) : dt_(dt) {

    }

    vector_t RKIntegrator::CalcIntegral(const vector_t& ic, const Trajectory& traj, double init_time,
                                        int num_steps,
                                        SingleRigidBodyModel& model, const vector_t& ref_state) {
        // TODO: DMA
        vector_t state = ic;
        for (int i = 0; i < num_steps; i++) {
            // TODO: DMA
            vector_t k1 = model.CalcDynamics(state, traj, init_time, ref_state);
            vector_t y1 = state + k1 * (dt_ / 2.0);
            vector_t k2 = model.CalcDynamics(y1, traj, init_time + dt_ / 2.0, ref_state);
            state = state + k2*dt_;
        }

        return state;
    }

    void RKIntegrator::CalcDerivWrtStateSingleStep(const matrix_t& dfdx, matrix_t& A) {
        // TODO: .noalias()
        A = matrix_t::Identity(dfdx.rows(), dfdx.cols()) + (dt_*dt_/2.0)*(dfdx*dfdx) + dt_*dfdx;
    }

    void RKIntegrator::CalcDerivWrtInputSingleStep(const matrix_t& dfdu, const matrix_t& dfdx,
                                                       matrix_t& B) {
        // TODO: .noalias()
        B = dt_*dfdu + (dt_*dt_/2.0)*(dfdx*dfdu);
    }

    void RKIntegrator::CalcLinearTermDiscretization(const vector_t& C, const vector_t& C2, const matrix_t& A,
                                                        vector_t& C_out) {
        // TODO: .noalias()
        C_out = dt_*C2 + (dt_*dt_/2.0)*A*C;
    }

    double RKIntegrator::GetDt() const {
        return dt_;
    }

    void RKIntegrator::DiscretizeLinearDynamics(matrix_t& A, matrix_t& B, vector_t& C, const vector_t& C2) {
        const matrix_t dfdx(A);    // TODO: DMA
        const matrix_t dfdu(B);    // TODO: DMA
        const vector_t constant(C);    // TODO: DMA

        CalcDerivWrtStateSingleStep(dfdx, A);
        CalcDerivWrtInputSingleStep(dfdu, dfdx, B);    // TODO: Technically this Ac should be evaluated at a different time
        CalcLinearTermDiscretization(constant, C2, dfdx, C);
    }
}