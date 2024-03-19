//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "integrator.h"

namespace mpc {
    Integrator::Integrator(double dt) : dt_(dt) {}

    double Integrator::GetDt() const {
        return dt_;
    }
}