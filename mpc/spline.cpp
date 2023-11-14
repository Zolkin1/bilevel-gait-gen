//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <cmath>
#include <cassert>

#include "include/spline.h"

namespace mpc {
    Spline::Spline(int num_polys, const std::vector<double> &times, bool start_on_poly) {
        // Determine the total number of polynomials. Depends on if we start on a single or a multiple
        if (start_on_poly) {
            total_poly_ = num_polys * (std::ceil((times.size() - 1)/2) + 1) + std::floor((times.size() - 1)/2) + 1;
        } else {
            total_poly_ = num_polys * (std::floor((times.size() - 1)/2) + 1) + std::ceil((times.size() - 1)/2) + 1;
        }
        for (int i = 0; i < total_poly_; i++) {
            poly_vars_.push_back(ZERO_POLY);
        }

        if (times.size() > 1) {
            for (int i = 1; i < times.size(); i++) {
                if (times.at(i - 1) >= times.at(i)) {
                    throw std::runtime_error("Switching times for the spline must be strictly increasing.");
                }
            }
        }
        poly_times_ = times;
    }

    double Spline::ValueAt(double time) const {
        assert(time >= 0);
        for (int i = 0; i < poly_times_.size(); i++) {
            if (time <= poly_times_.at(i)) {
                double DeltaT = -1;
                if (i == 0) {
                    DeltaT = poly_times_.at(i);
                } else {
                    DeltaT = poly_times_.at(i) - poly_times_.at(i-1);
                }
                return EvalPoly(poly_vars_.at(i), time - poly_times_.at(i), DeltaT);
            }
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline.");
    }

    void Spline::SetPolyVars(int poly_num, const std::array<double, POLY_ORDER>& vars) {
        poly_vars_.at(poly_num) = vars;
    }

    double Spline::EvalPoly(const std::array<double, POLY_ORDER>& poly_vals, double time, double DeltaT) {
        double a2 = -(1/pow(DeltaT, 2))*3*(poly_vals.at(0) - poly_vals.at(1)) -
                (1/DeltaT)*(2*poly_vals.at(2) + poly_vals.at(3));
        double a3 = (1/pow(DeltaT, 3))*2*(poly_vals.at(0) - poly_vals.at(1)) +
                (1/pow(DeltaT, 2))*(poly_vals.at(2) + poly_vals.at(3));
        double val = poly_vals.at(0) + poly_vals.at(1)*time + a2*pow(time, 2) + a3*pow(time, 3);

        return val;
    }

    int Spline::GetTotalPoly() const {
        return total_poly_;
    }

} // mpc