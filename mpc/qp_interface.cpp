//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "include/qp_interface.h"

namespace mpc {

    QPInterface::QPInterface() {}

    vector_t QPInterface::GetInfinity(int size) const {
        vector_t infty = vector_t::Constant(size, 1e30);

        return infty;
    }

} // mpc