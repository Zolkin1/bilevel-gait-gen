//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
#define BILEVEL_GAIT_GEN_OSQP_INTERFACE_H

#include <Eigen/SparseCore>

#include "OsqpEigen/OsqpEigen.h"

#include "qp_interface.h"
#include "timer.h"

namespace mpc {
    class OSQPInterface : public QPInterface {
    public:
        OSQPInterface(QPData data, bool verbose);

        void SetupQP(QPData& data, const vector_t& warm_start) override;

        // TODO: remove data
        vector_t Solve(const QPData& data) override;

        vector_t GetInfinity(int size) const override;

        std::string GetSolveQuality() const override;

        vector_t GetDualSolution() const override;

        void ConfigureForInitialRun() const override;

        void ConfigureForRealTime(double run_time_iters) const override;



    protected:
    private:
//        void ConvertDataToOSQPConstraints(const QPData& data);

//        void ConvertDataToOSQPCost(const QPData& data);

        // QP Solver
        OsqpEigen::Solver qp_solver_;

        vector_t prev_dual_sol_;

        int run;

        bool verbose_;

        utils::Timer osqp_interface_timer_;
        utils::Timer sparse_conversion_timer_;
        utils::Timer matrix_gather_timer_;

        int prev_num_constraints_;
        int prev_num_decision_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_OSQP_INTERFACE_H
