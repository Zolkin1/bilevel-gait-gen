//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include "qp/qp_interface.h"

namespace mpc {

    QPInterface::QPInterface(int num_decision_vars) {
        prev_qp_sol_ = vector_t::Zero(num_decision_vars);
    }

    vector_t QPInterface::GetInfinity(int size) const {
        vector_t infty = vector_t::Constant(size, 1e30);

        return infty;
    }

    std::string QPInterface::GetSolveQualityAsString() const {
        switch (GetSolveQuality()) {
            case Solved:
                return "Solved";
            case SolvedInacc:
                return "Solved Inaccurate";
            case PrimalInfeasible:
                return "Primal Infeasible";
            case DualInfeasible:
                return "Dual Infeasible";
            case PrimalInfeasibleInacc:
                return "Primal Infeasible Inaccurate";
            case DualInfeasibleInacc:
                return "Dual Infeasible Inaccurate";
            case MaxIter:
                return "Max Iter Reached";
            case Unsolved:
                return "Unsolved";
            default:
                return "Other";
        }
    }

} // mpc