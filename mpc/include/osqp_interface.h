//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
#define BILEVEL_GAIT_GEN_OSQP_INTERFACE_H

#include "OsqpEigen/OsqpEigen.h"

#include "qp_interface.h"
#include "timer.h"

namespace mpc {
    class OSQPInterface : public QPInterface {
    public:
        OSQPInterface(QPData data, bool verbose);

        void SetupQP(const QPData& data, const vector_t& warm_start) override;

        // TODO: remove data
        vector_t Solve(const QPData& data) override;

        vector_t GetInfinity(int size) const override;

        std::string GetSolveQuality() const override;

        vector_t GetDualSolution() const override;

        void ConfigureForInitialRun() const override;

        void ConfigureForRealTime() const override;

    protected:
    private:
        void ConvertDataToOSQPConstraints(const QPData& data);

        void ConvertDataToOSQPCost(const QPData& data);

        // QP Solver
        OsqpEigen::Solver qp_solver_;

        // OSQP specific vectors and matricies
        matrix_t A_, P_;
//        Eigen::SparseMatrix<double> A_, P_;
        vector_t lb_, ub_, w_;

        vector_t prev_dual_sol_;

        int run;

        bool verbose_;

        utils::Timer osqp_interface_timer_;
        utils::Timer sparse_conversion_timer_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
