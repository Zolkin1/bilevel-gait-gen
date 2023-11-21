//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//
#include <iostream>

#include "euler_integrator.h"
#include "centroidal_model.h"

namespace mpc {
    EulerIntegrator::EulerIntegrator(double dt) : Integrator(dt) {}

    vector_t EulerIntegrator::CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time,
                                           double final_time, const CentroidalModel& model) {
        if (final_time < init_time) {
            throw std::runtime_error("Final time for integration less than initial time.");
        }

        if ((final_time - init_time)/dt_ != round((final_time - init_time)/dt_)) {
            std::cerr << "Integration time not evenly divisible by dt. Rounding to the nearest divisible time." << std::endl;
        }

        int num_iters = static_cast<int>(round((final_time - init_time)/dt_));
        vector_t val = ic;
        double time = init_time;
        for (int i = 0; i < num_iters; i++) {
            val += model.CalcDynamics(val, input, time)*dt_;
        }

        return val;
    }

    vector_t EulerIntegrator::CalcIntegral(const mpc::vector_t& ic, const mpc::Inputs& input, double init_time,
                                           int num_steps, const mpc::CentroidalModel& model) {
        vector_t val = ic;
        double time = init_time;
        for (int i = 0; i < num_steps; i++) {
            val += model.CalcDynamics(val, input, time)*dt_;
        }

        return val;
    }

    matrix_t EulerIntegrator::CalcDerivWrtStateSingleStep(const mpc::vector_t &ic, const mpc::matrix_t &dfdx) {
        return dfdx*dt_ + matrix_t::Identity(ic.size() - 1, ic.size() - 1);
    }

    matrix_t EulerIntegrator::CalcDerivWrtInputSingleStep(const mpc::vector_t &ic, const mpc::matrix_t &dfdu) {
        return dfdu * dt_;
    }
}