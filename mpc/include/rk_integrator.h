//
// Created by zolkin on 12/8/23.
//

#ifndef BILEVEL_GAIT_GEN_RK_INTEGRATOR_H
#define BILEVEL_GAIT_GEN_RK_INTEGRATOR_H

#include "integrator.h"
//#include "centroidal_model.h"

namespace mpc {
    class RKIntegrator {
    public:
        RKIntegrator(double dt);

        vector_t CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time, double final_time,
                              CentroidalModel& model);

        vector_t CalcIntegral(const mpc::vector_t &ic, const mpc::Inputs &input, double init_time, int num_steps,
                              CentroidalModel& model, const vector_t& ref_state);

        // Only does a single step
        void CalcDerivWrtStateSingleStep(const vector_t& ic, const matrix_t& dfdx,
                                             matrix_t& A);

        // Only does a single step
        void CalcDerivWrtInputSingleStep(const vector_t& ic, const matrix_t& dfdu,
                                             const matrix_t& dfdx,
                                             matrix_t& B);

        void CalcLinearTermDiscretization(const vector_t& C, const vector_t& C2, const matrix_t& A,
                                          vector_t& C_out);

        double GetDt() const;
    protected:
    private:
        double dt_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_RK_INTEGRATOR_H
