//
// Created by zolkin on 12/8/23.
//

#include "rk_integrator.h"
#include "centroidal_model.h"

namespace mpc {
    RKIntegrator::RKIntegrator(double dt) : dt_(dt) {

    }

    vector_t RKIntegrator::CalcIntegral(const vector_t& ic, const Inputs& input, double init_time, int num_steps,
                                const CentroidalModel& model, const vector_t& ref_state) {
        vector_t state = ic;
        for (int i = 0; i < num_steps; i++) {
            vector_t k1 = model.CalcDynamics(state, input, init_time, ref_state);
            vector_t y1 = state + k1 * (dt_ / 2.0);
            vector_t k2 = model.CalcDynamics(y1, input, init_time + dt_ / 2.0, ref_state);
            state = state + k2*dt_;
        }

        return state;
    }

    matrix_t RKIntegrator::CalcDerivWrtStateSingleStep(const vector_t& ic, const matrix_t& dfdx) {
        return matrix_t::Identity(dfdx.rows(), dfdx.cols()) + (dt_*dt_/2)*(dfdx*dfdx) + dt_*dfdx;
    }

    matrix_t RKIntegrator::CalcDerivWrtInputSingleStep(const vector_t& ic, const matrix_t& dfdu, const matrix_t& dfdx) {
        return (dt_*dt_/2)*(dfdx*dfdu) + dt_*dfdu;
    }

    vector_t RKIntegrator::CalcLinearTermDiscretization(const vector_t& C, const vector_t& C2, const matrix_t& A) {
        return dt_*C2 + (dt_*dt_/2)*A*C;
    }

    double RKIntegrator::GetDt() const {
        return dt_;
    }
}