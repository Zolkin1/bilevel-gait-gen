//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//
#include <iostream>

#include "euler_integrator.h"
#include "models/centroidal_model.h"

namespace mpc {
    EulerIntegrator::EulerIntegrator(double dt) : Integrator(dt) {}

    vector_t EulerIntegrator::CalcIntegral(const mpc::vector_t &ic, const mpc::Trajectory& traj,
                                           double init_time, double final_time, SingleRigidBodyModel& model) {
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
            val += model.CalcDynamics(val, traj, time, ic)*dt_;
        }

        return val;
    }

    vector_t EulerIntegrator::CalcIntegral(const vector_t& ic, const Trajectory& traj, double init_time,
                                           int num_steps, SingleRigidBodyModel& model, const vector_t& ref_state) {
        vector_t val = ic;
        double time = init_time;
        for (int i = 0; i < num_steps; i++) {
            val += model.CalcDynamics(val, traj, time, ref_state)*dt_;
        }

        return val;
    }

    void EulerIntegrator::CalcDerivWrtStateSingleStep(const mpc::vector_t &ic, const mpc::matrix_t &dfdx,
                                                          Eigen::Ref<matrix_t> A) {
        A.noalias() = dfdx*dt_ + matrix_t::Identity(ic.size() - 1, ic.size() - 1);
    }

    void EulerIntegrator::CalcDerivWrtInputSingleStep(const mpc::vector_t &ic, const mpc::matrix_t &dfdu,
                                                          const matrix_t& dfdx, Eigen::Ref<matrix_t> B) {
        B.noalias() = dfdu * dt_;
    }
}