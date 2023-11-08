//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_QP_INTERFACE_H
#define BILEVEL_GAIT_GEN_QP_INTERFACE_H

#include "mpc.h"

namespace mpc {
    /**
     * Base class to handle all interfaces to the underlying QP solver.
     * Each QP solver should have its own child class that implements these functions.
     */
    class QPInterface {
    public:
        QPInterface();

        void SetupQP();

        Trajectory SolveQP;
    protected:
    private:


    };
}


#endif //BILEVEL_GAIT_GEN_QP_INTERFACE_H
