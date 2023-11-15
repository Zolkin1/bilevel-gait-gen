//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
#define BILEVEL_GAIT_GEN_OSQP_INTERFACE_H

#include "OsqpEigen/OsqpEigen.h"

#include "qp_interface.h"


namespace mpc {
    class OSQPInterface : public QPInterface {
    public:
        OSQPInterface();

        void SetupQP(const QPData& data) override;

        vector_t Solve() override;

        vector_t GetInfinity(int size) const override;
    protected:
    private:
        // TODO: write
        void ConvertDataToOSQPConstraints(const QPData& data);

        void ConvertDataToOSQPCost(const QPData& data);

        // QP Solver
        OsqpEigen::Solver qp_solver_;

        // OSQP specific vectors and matricies
        matrix_t A_, P_;
        vector_t lb_, ub_, w_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
