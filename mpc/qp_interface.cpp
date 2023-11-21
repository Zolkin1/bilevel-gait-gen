//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "include/qp_interface.h"

namespace mpc {

    QPInterface::QPInterface(int num_decision_vars) {
        prev_qp_sol_ = vector_t::Zero(num_decision_vars);
    }

    vector_t QPInterface::GetInfinity(int size) const {
        vector_t infty = vector_t::Constant(size, 1e30);

        return infty;
    }

} // mpc