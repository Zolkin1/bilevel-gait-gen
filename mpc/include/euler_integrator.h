//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_EULER_INTEGRATOR_H
#define BILEVEL_GAIT_GEN_EULER_INTEGRATOR_H

#include "integrator.h"

namespace mpc {
    class EulerIntegrator : public Integrator {
    public:
        EulerIntegrator(double dt);

        vector_t CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time, double final_time,
                              CentroidalModel& model) override;

        vector_t CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time, int num_steps,
                              CentroidalModel& model, const vector_t& ref_state) override;

        // Only does a single step
        void CalcDerivWrtStateSingleStep(const vector_t& ic, const matrix_t& dfdx,
                                         Eigen::Ref<matrix_t> A) override;

        // Only does a single step
        void CalcDerivWrtInputSingleStep(const vector_t& ic, const matrix_t& dfdu,
                                             const matrix_t& dfdx, Eigen::Ref<matrix_t> B) override;
    protected:
    private:
    };
}


#endif //BILEVEL_GAIT_GEN_EULER_INTEGRATOR_H
