//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <cmath>
#include <cassert>

#include "spline.h"

namespace mpc {
    Spline::Spline(int num_polys, const std::vector<double> &times, bool start_on_poly) :
                    ZERO_CONSTANT(1,0), ZERO_POLY(2,0), start_on_constant_(!start_on_poly), num_polys_(num_polys) {

        poly_times_.push_back(0);

        for (int time_idx = 0; time_idx < times.size(); time_idx++) {
            if (start_on_constant_) {
                if ((time_idx) % 2 == 0) {      // On a constant
                    poly_vars_.push_back(ZERO_CONSTANT);
                    poly_vars_.push_back(ZERO_CONSTANT);
                    num_constant_++;
                    total_poly_++;
                    num_all_poly_vars++;

                    poly_times_.push_back(times.at(time_idx));
                } else {        // On a polynomial
                    for (int i = 0; i < num_polys_ - 1; i++) {
                        poly_vars_.push_back(ZERO_POLY);

                        poly_times_.push_back((times.at(time_idx) - times.at(time_idx - 1))/num_polys
                        + times.at(time_idx - 1));

                        num_all_poly_vars += 2;
                    }

                }
            } else {
                if ((time_idx + 1) % 2 == 0) {      // On a constant
                    poly_vars_.push_back(ZERO_CONSTANT);
                    poly_vars_.push_back(ZERO_CONSTANT);
                    num_constant_++;
                    total_poly_++;
                    num_all_poly_vars++;

                    poly_times_.push_back(times.at(time_idx));
                } else {        // On a polynomial
                    total_poly_ += num_polys_;
                    for (int i = 0; i < num_polys_ - 1; i++) {
                        poly_vars_.push_back(ZERO_POLY);

                        if (time_idx == 0) {
                            poly_times_.push_back(times.at(time_idx)/num_polys);
                        } else {
                            poly_times_.push_back((times.at(time_idx) - times.at(time_idx - 1))/num_polys
                                                  + times.at(time_idx - 1));
                        }
                        num_all_poly_vars += 2;
                    }
                }
            }
        }
    }

    double Spline::ValueAt(double time) const {
        std::vector<double> vars = GetPolyVars(time);

        int i = GetPolyIdx(time);
        double DeltaT = poly_times_.at(i) - poly_times_.at(i-1);
        time = time - poly_times_.at(i-1);

        return EvalPoly(vars, time, DeltaT);
    }

    // TODO: Change to time based! Does not make sense as written!
    void Spline::SetPolyVars(int poly_num, const std::vector<double>& vars) {
        assert(poly_num > 0);

        if (vars.size() != poly_vars_.at(poly_num).size()) {
            throw std::runtime_error("Provided polynomial does not match the order of the current polynomial.");
        }

        poly_vars_.at(poly_num) = vars;

        // Check for constants
        if (poly_vars_.at(poly_num).size() == 1 && poly_vars_.at(poly_num-1).size() == 1) {
            poly_vars_.at(poly_num-1) = vars;
        } else if (poly_vars_.at(poly_num).size() == 1 && poly_vars_.at(poly_num+1).size() == 1) {
            poly_vars_.at(poly_num+1) = vars;
        }
    }

    void Spline::SetAllSplineVars(const std::vector<std::vector<double>>& vars) {
        if (vars.size() != poly_vars_.size()) {
            throw std::runtime_error("Not enough polynomial sections provided.");
        }

        for (int i = 0; i < vars.size(); i++) {
            if (vars.at(i).size() != poly_vars_.at(i).size()) {
                throw std::runtime_error("One of the provided variable vectors does not have the correct size.");
            }
        }

        poly_vars_ = vars;
    }

    // Time between 0 and DT
    double Spline::EvalPoly(const std::vector<double>& vars, double time, double DeltaT) {
        double a2 = -(1 / pow(DeltaT, 2)) * 3 * (vars.at(0) - vars.at(1)) -
                    (1 / DeltaT) * (2 * vars.at(2) + vars.at(3));
        double a3 = (1 / pow(DeltaT, 3)) * 2 * (vars.at(0) - vars.at(1)) +
                    (1 / pow(DeltaT, 2)) * (vars.at(2) + vars.at(3));
        double val = vars.at(0) + vars.at(1) * time + a2 * pow(time, 2) + a3 * pow(time, 3);

        return val;
    }

    int Spline::GetTotalPolyVars() const {
        return num_all_poly_vars;
    }

    int Spline::GetTotalPoly() const {
        return total_poly_;
    }

    std::vector<double> Spline::GetPolyVars(double time) const {
        int i = GetPolyIdx(time);

        if (IsConstantPoly(i)) {
            // Constant
            std::vector<double> vars{poly_vars_.at(i).at(0)};
            return vars;
        }

        // Polynomial; gather the variables that describe it
        std::vector<double> vars(4,0);
        vars.at(0) = poly_vars_.at(i-1).at(0);      // x0
        vars.at(1) = poly_vars_.at(i).at(0);        // x1

        // Check if the "end" point is a constant, x1dot
        if (poly_vars_.at(i).size() == 1) {
            vars.at(4) = 0;
        } else {
            vars.at(4) = poly_vars_.at(i).at(1);
        }

        // Check if the "start" point is a constant, x0dot
        if (poly_vars_.at(i-1).size() == 1) {
            vars.at(3) = 0;
        } else {
            vars.at(3) = poly_vars_.at(i).at(1);
        }

        return vars;
    }

    vector_t Spline::GetAllPolyVars() const {
        vector_t all_vars = vector_t::Zero(num_all_poly_vars);
        int idx = 0;
        for (int poly = 0; poly < poly_vars_.size(); poly++) {
            for (int j = 0; j < poly_vars_.at(poly).size(); j++) {
                all_vars(idx + j) = poly_vars_.at(poly).at(j);
            }
            idx += poly_vars_.at(poly).size();

            // Constants are always followed by redundant constants
            if (poly_vars_.at(poly).size() == 1) {
                poly++;
            }
        }

        return all_vars;
    }

    int Spline::GetPolyIdx(double time) const {
        assert(time >= 0);
        for (int i = 1; i < poly_times_.size(); i++) {     // Start at i = 1 because t = 0  at 0
            if (time <= poly_times_.at(i)) {
                return i;
            }
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline. Time too large.");
    }

    vector_t Spline::GetPolyVarsLin(double time) const {
        int i = GetPolyIdx(time);

        if (IsConstantPoly(i)) {
            vector_t vars(1);
            vars << 1;
            return vars;
        }

        // Polynomial
        double DeltaT = poly_times_.at(i) - poly_times_.at(i-1);
        time = time - poly_times_.at(i-1);

        if (poly_vars_.at(i).size() == 1) {
            // Final time derivative is fixed at 0
            vector_t vars = vector_t::Zero(3);
            vars(0) = 1 - 1/pow(DeltaT, 2)*3*pow(time,2) + 1/pow(DeltaT, 3)*2*pow(time, 3);     // x0 coef
            vars(1) = 1/pow(DeltaT, 2)*3*pow(time,2) - 1/pow(DeltaT, 3)*2*pow(time, 3);         // x1 coef

            vars(2) = time - 1/DeltaT*2*pow(time, 2) + 1/pow(DeltaT, 2)*pow(time, 3);           // x0dot coef
            return vars;
        } else if (poly_vars_.at(i-1).size() == 1) {
            // Initial time derivative is fixed at 0
            vector_t vars = vector_t::Zero(3);
            vars(0) = 1 - 1/pow(DeltaT, 2)*3*pow(time,2) + 1/pow(DeltaT, 3)*2*pow(time, 3);     // x0 coef
            vars(1) = 1/pow(DeltaT, 2)*3*pow(time,2) - 1/pow(DeltaT, 3)*2*pow(time, 3);         // x1 coef

            vars(2) = -1/DeltaT*pow(time, 2) + 1/pow(DeltaT, 2)*pow(time, 3);                   // x1dot coef

            return vars;
        } else {
            vector_t vars = vector_t::Zero(4);
            vars(0) = 1 - 1/pow(DeltaT, 2)*3*pow(time,2) + 1/pow(DeltaT, 3)*2*pow(time, 3);     // x0 coef
            vars(1) = 1/pow(DeltaT, 2)*3*pow(time,2) - 1/pow(DeltaT, 3)*2*pow(time, 3);         // x1 coef

            vars(2) = time - 1/DeltaT*2*pow(time, 2) + 1/pow(DeltaT, 2)*pow(time, 3);           // x0dot coef
            vars(3) = -1/DeltaT*pow(time, 2) + 1/pow(DeltaT, 2)*pow(time, 3);                   // x1dot coef
            return vars;
        }

    }

    bool Spline::IsConstantPoly(int idx) const {
        return (poly_vars_.at(idx).size() == 1 && poly_vars_.at(idx-1).size() == 1);
    }

} // mpc