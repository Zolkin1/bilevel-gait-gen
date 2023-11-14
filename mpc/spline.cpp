//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <cmath>
#include <cassert>

#include "spline.h"

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
                    time = time - poly_times_.at(i-1);
                }
                return EvalPoly(poly_vars_.at(i), time, DeltaT);
            }
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline. Time too large.");
    }

    void Spline::SetPolyVars(int poly_num, const std::array<double, POLY_ORDER>& vars) {
        poly_vars_.at(poly_num) = vars;
    }

    // Time between 0 and 1
    double Spline::EvalPoly(const std::array<double, POLY_ORDER>& poly_vals, double time, double DeltaT) {
        double a2 = -(1/pow(DeltaT, 2))*3*(poly_vals.at(0) - poly_vals.at(1)) -
                (1/DeltaT)*(2*poly_vals.at(2) + poly_vals.at(3));
        double a3 = (1/pow(DeltaT, 3))*2*(poly_vals.at(0) - poly_vals.at(1)) +
                (1/pow(DeltaT, 2))*(poly_vals.at(2) + poly_vals.at(3));
        double val = poly_vals.at(0) + poly_vals.at(1)*time + a2*pow(time, 2) + a3*pow(time, 3);

        return val;
    }

    int Spline::GetTotalPolyVars() const {
        return total_poly_*POLY_ORDER;
    }

    int Spline::GetTotalPoly() const {
        return total_poly_;
    }

    Eigen::Vector<double, Spline::POLY_ORDER> Spline::GetPolyVars(double time) const {
        assert(time >= 0);
        for (int i = 0; i < poly_times_.size(); i++) {
            if (time <= poly_times_.at(i)) {
                Eigen::Vector<double, POLY_ORDER> vars = Eigen::Vector<double, POLY_ORDER>::Zero();
                for (int j = 0; j < poly_vars_.at(i).size(); j++) {
                    vars(j) = poly_vars_.at(i).at(j);
                }
                return vars;
            }
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline. Time too large.");
    }

    vector_t Spline::GetAllPolyVars() const {
        vector_t all_vars = vector_t::Zero(poly_vars_.size()*POLY_ORDER);
        for (int i = 0; i < poly_vars_.size(); i++) {
            for (int j = 0; j < POLY_ORDER; j++) {
                all_vars(i*POLY_ORDER + j) = poly_vars_.at(i).at(j);
            }
        }

        return all_vars;
    }

    int Spline::GetPolyIdx(double time) const {
        assert(time >= 0);
        for (int i = 0; i < poly_times_.size(); i++) {
            if (time <= poly_times_.at(i)) {
                return i;
            }
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline. Time too large.");
    }

    Eigen::Vector<double, Spline::POLY_ORDER> Spline::GetPolyVarsLin(double time) const {
        assert(time >= 0);
        for (int i = 0; i < poly_times_.size(); i++) {
            if (time <= poly_times_.at(i)) {
                double DeltaT = -1;
                if (i == 0) {
                    DeltaT = poly_times_.at(i);
                } else {
                    DeltaT = poly_times_.at(i) - poly_times_.at(i-1);
                    time = time - poly_times_.at(i-1);
                }
                Eigen::Vector<double, POLY_ORDER> vars = Eigen::Vector<double, POLY_ORDER>::Zero();
                vars(0) = 1 - 1/pow(DeltaT, 2)*3*pow(time,2) + 1/pow(DeltaT, 3)*2*pow(time, 3);
                vars(1) = time - 1/DeltaT*2*pow(time, 2) + 1/pow(DeltaT, 2)*pow(time, 3);
                vars(2) = 1/pow(DeltaT, 2)*3*pow(time,2) - 1/pow(DeltaT, 3)*2*pow(time, 3);
                vars(3) = -1/DeltaT*pow(time, 2) + 1/pow(DeltaT, 2)*pow(time, 3);
                return vars;
            }
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline. Time too large.");
    }

} // mpc