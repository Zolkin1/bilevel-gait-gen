//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_INTEGRATOR_H
#define BILEVEL_GAIT_GEN_INTEGRATOR_H

#include <Eigen/Core>

#include "inputs.h"

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using matrix_t = Eigen::MatrixXd;

    class CentroidalModel;

    class Integrator {
    public:
        Integrator(double dt);

        virtual vector_t CalcIntegral(const vector_t& ic, const Inputs& input, double init_time, double final_time,
                                      const CentroidalModel& model) = 0;

        virtual matrix_t CalcDerivWrtStateSingleStep(const vector_t& ic, const matrix_t& dfdx) = 0;

        virtual matrix_t CalcDerivWrtInputSingleStep(const vector_t& ic, const matrix_t& dfdu) = 0;

    protected:
        double dt_;

    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_INTEGRATOR_H
