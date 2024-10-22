//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_QP_INTERFACE_H
#define BILEVEL_GAIT_GEN_QP_INTERFACE_H

#include "trajectory.h"
#include "qp_data.h"

namespace mpc {
    enum SolveQuality {
        Solved = 0,
        SolvedInacc = 1,
        MaxIter = 2,
        PrimalInfeasible = 3,
        DualInfeasible = 4,
        PrimalInfeasibleInacc = 5,
        DualInfeasibleInacc = 6,
        Unsolved = 7,
        Other = 8
    };

    using vector_t = Eigen::VectorXd;

    /**
     * Base class to handle all interfaces to the underlying QP solver.
     * Each QP solver should have its own child class that implements these functions.
     */
    class QPInterface {
    public:

        QPInterface(int num_decision_vars);

        /*
         * Notes:
         * I tried to get the Intel MKL library to work so I could try that solver. It did not seem to work,
         * it always claimed it couldn't find the openmp library, which might be because it was trying to use intel
         * openmp but I am using gnu openmp. Might be worth coming back to this to see if the multi-threaded speed up is
         * large.
         *
         * OSQP eigen never seems to let me use a time limit because it doesn't set the profiling variable even though
         * OSQP is built with profiling. Can probably fix this by either (1) opening and issue on their github, (2)
         * modifying source, or (3) finding another work around. Might want to use this to ensure timing requirements.
         */
        virtual void SetupQP(QPData& data, const vector_t& warm_start) = 0;

        virtual vector_t Solve(const QPData& data) = 0;

        virtual vector_t GetInfinity(int size) const;

        virtual SolveQuality GetSolveQuality() const = 0;

        std::string GetSolveQualityAsString() const;

        virtual vector_t GetDualSolution() const = 0;

        virtual void ConfigureForInitialRun() = 0;

        virtual void ConfigureForRealTime(double run_time_iters) = 0;

    protected:
        vector_t prev_qp_sol_;
    private:
    };
}


#endif //BILEVEL_GAIT_GEN_QP_INTERFACE_H
