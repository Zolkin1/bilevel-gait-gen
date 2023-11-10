//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_EULER_INTEGRATOR_H
#define BILEVEL_GAIT_GEN_EULER_INTEGRATOR_H

#include "integrator.h"

namespace mpc {
    class EulerIntegrator : public Integrator {
        EulerIntegrator(double dt);

        vector_t CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time, double final_time,
                              const CentroidalModel& model) override;

        // Only does a single step
        matrix_t CalcDerivWrtStateSingleStep(const vector_t& ic, const matrix_t& dfdx) override;

        // Only does a single step
        matrix_t CalcDerivWrtInputSingleStep(const vector_t& ic, const matrix_t& dfdu) override;
    };
}


#endif //BILEVEL_GAIT_GEN_EULER_INTEGRATOR_H
