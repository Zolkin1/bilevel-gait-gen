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

        Trajectory Solve() override;
    protected:
    private:
        // TODO: write
        matrix_t ConvertDataToOSQPConstraintMat(const QPData& data);
        vector_t ConvertDataToOSQPlb(const QPData& data);
        vector_t ConvertDataToOSQPub(const QPData& data);

        matrix_t ConvertDataToOSQPCostMat(const QPData& data);
        matrix_t ConvertDataToOSQPCostVec(const QPData& data);

        // QP Solver
        OsqpEigen::Solver qp_solver_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
