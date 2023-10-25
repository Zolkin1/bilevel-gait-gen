//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "controller.h"

namespace simulator {
    Controller::Controller(double control_rate) : rate_(1.0/control_rate) {}

    double Controller::GetRate() const {
        return rate_;
    }
} // simulator