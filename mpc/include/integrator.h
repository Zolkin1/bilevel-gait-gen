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

        double GetDt() const;
        virtual vector_t CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time, double final_time,
                              CentroidalModel& model) = 0;

        virtual vector_t CalcIntegral(const vector_t& ic, const Inputs& input, double init_time, int num_steps,
                                      CentroidalModel& model, const vector_t& ref_state) = 0;

        virtual void CalcDerivWrtStateSingleStep(const vector_t& ic, const matrix_t& dfdx,
                                                 Eigen::Ref<matrix_t> A) = 0;

        virtual void CalcDerivWrtInputSingleStep(const vector_t& ic, const matrix_t& dfdu,
                                                     const matrix_t& dfdx,
                                                     Eigen::Ref<matrix_t> B) = 0;
    protected:
        double dt_;

    private:
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_INTEGRATOR_H
