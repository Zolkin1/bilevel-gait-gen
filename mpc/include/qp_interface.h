//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_QP_INTERFACE_H
#define BILEVEL_GAIT_GEN_QP_INTERFACE_H

#include "trajectory.h"
#include "qp_data.h"

namespace mpc {

    using vector_t = Eigen::VectorXd;

    /**
     * Base class to handle all interfaces to the underlying QP solver.
     * Each QP solver should have its own child class that implements these functions.
     */
    class QPInterface {
    public:
        QPInterface();

        virtual void SetupQP(const QPData& data) = 0;

        virtual Trajectory Solve() = 0;
    protected:
    private:


    };
}


#endif //BILEVEL_GAIT_GEN_QP_INTERFACE_H
