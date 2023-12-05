//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <cmath>
#include <cassert>

#include "spline.h"

namespace mpc {
    Spline::Spline(int num_polys, const std::vector<double> &times, bool start_on_poly) :
                    start_on_constant_(!start_on_poly), num_polys_(num_polys) {

        std::vector<double> const ZERO_CONSTANT = {0};
        std::vector<double> const ZERO_POLY = {0, 0};

        poly_times_.push_back(0);
        num_all_poly_vars = 0;
        total_poly_ = 0;
        num_constant_ = 0;

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
                    total_poly_ += num_polys_;

                    for (int i = 0; i < num_polys_ - 1; i++) {
                        poly_vars_.push_back(ZERO_POLY);

                        poly_times_.push_back((times.at(time_idx) - times.at(time_idx - 1))/num_polys
                        + poly_times_.at(poly_times_.size()-1));

                        num_all_poly_vars += 2;
                    }

                    poly_times_.push_back(times.at(time_idx));
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
                    if (time_idx == 0) {
                        poly_times_.push_back(times.at(time_idx)/num_polys);
                        poly_vars_.push_back(ZERO_CONSTANT);
                        num_all_poly_vars++;
                        num_constant_++;
                    }

                    for (int i = 0; i < num_polys_ - 1; i++) {
                        if (time_idx == 0) {
                            poly_times_.push_back(times.at(time_idx)/num_polys + poly_times_.at(poly_times_.size()-1));
                        } else {
                            poly_times_.push_back((times.at(time_idx) - times.at(time_idx - 1)) / num_polys
                                                  + poly_times_.at(poly_times_.size()-1));
                        }

                        poly_vars_.push_back(ZERO_POLY);
                        num_all_poly_vars += 2;
                    }

//                    poly_times_.push_back(times.at(time_idx));
                }
            }
        }

        // Add final elements
        if ((start_on_constant_ && (times.size()) % 2 == 0) || (!start_on_constant_ && (times.size()-1) % 2 == 0)) {
//            poly_times_.push_back(times.at(times.size()-1));
            poly_vars_.push_back(ZERO_CONSTANT);
            num_all_poly_vars += 1;
            num_constant_++;
        } else {
//            poly_vars_.push_back(ZERO_POLY);
//            num_all_poly_vars += 2;
        }
    }

    double Spline::ValueAt(double time) const {
        std::array<double, POLY_ORDER> vars = GetPolyVars(time);

        int i = GetPolyIdx(time);
        double DeltaT = poly_times_.at(i) - poly_times_.at(i-1);
        time = time - poly_times_.at(i-1);

        return EvalPoly(vars, time, DeltaT);
    }

    void Spline::SetPolyVars(int poly_time, const std::vector<double>& vars) {
        assert(poly_time >= 0);
        assert(poly_time <= poly_times_.size());

        if (vars.size() != poly_vars_.at(poly_time).size()) {
            throw std::runtime_error("Provided polynomial does not match the order of the current polynomial.");
        }

        poly_vars_.at(poly_time) = vars;

        if (poly_time == 0) {
            // Don't need to check prev val
            if (poly_vars_.at(poly_time).size() == 1 && poly_vars_.at(poly_time + 1).size() == 1) {
                poly_vars_.at(poly_time + 1) = vars;
            }
        } else if (poly_time == poly_times_.size()) {
            // Don't need to check after val
            if (poly_vars_.at(poly_time).size() == 1 && poly_vars_.at(poly_time - 1).size() == 1) {
                poly_vars_.at(poly_time - 1) = vars;
            }
        } else {
            // Check for constants
            if (poly_vars_.at(poly_time).size() == 1 && poly_vars_.at(poly_time - 1).size() == 1) {
                poly_vars_.at(poly_time - 1) = vars;
            } else if (poly_time+1<poly_vars_.size() && poly_vars_.at(poly_time).size() == 1 && poly_vars_.at(poly_time + 1).size() == 1) {
                poly_vars_.at(poly_time + 1) = vars;
            }
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
    double Spline::EvalPoly(const std::array<double, POLY_ORDER>& vars, double time, double DeltaT) {
        double a2 = -(1 / pow(DeltaT, 2)) * 3 * (vars.at(0) - vars.at(1)) -
                    (1 / DeltaT) * (2 * vars.at(2) + vars.at(3));
        double a3 = (1 / pow(DeltaT, 3)) * 2 * (vars.at(0) - vars.at(1)) +
                    (1 / pow(DeltaT, 2)) * (vars.at(2) + vars.at(3));

        double val = vars.at(0) + vars.at(2) * time + a2 * pow(time, 2) + a3 * pow(time, 3);

        return val;
    }

    int Spline::GetTotalPolyVars() const {
        return num_all_poly_vars;
    }

    int Spline::GetTotalPoly() const {
        return total_poly_;
    }

    std::array<double, Spline::POLY_ORDER> Spline::GetPolyVars(double time) const {
        int i = GetPolyIdx(time);

//        if (IsConstantPoly(i)) {
//            // Constant
//            std::vector<double> vars{poly_vars_.at(i).at(0)};
//            return vars;
//        }

        // Polynomial; gather the variables that describe it
        std::array<double, Spline::POLY_ORDER> vars = {0,0,0,0};
        vars.at(0) = poly_vars_.at(i-1).at(0);      // x0
        vars.at(1) = poly_vars_.at(i).at(0);        // x1

        // Check if the "end" point is a constant, x1dot
        if (poly_vars_.at(i).size() == 1) {
            vars.at(3) = 0;
        } else {
            vars.at(3) = poly_vars_.at(i).at(1);
        }

        // Check if the "start" point is a constant, x0dot
        if (poly_vars_.at(i-1).size() == 1) {
            vars.at(2) = 0;
        } else {
            vars.at(2) = poly_vars_.at(i-1).at(1);
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
            if (poly_vars_.at(poly).size() == 1 && poly < poly_vars_.size()-1 && poly_vars_.at(poly+1).size() == 1) {
                poly++;
            }
        }

        return all_vars;
    }

    int Spline::GetPolyIdx(double time) const {
        assert(time >= poly_times_.at(0));
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

        if (IsConstantPoly(i) || time == poly_times_.at(0)) {
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
            vars(0) = 1 - (1/pow(DeltaT, 2))*3*pow(time,2) + (1/pow(DeltaT, 3))*2*pow(time, 3);     // x0 coef
            vars(2) = (1/pow(DeltaT, 2))*3*pow(time,2) - (1/pow(DeltaT, 3))*2*pow(time, 3);         // x1 coef

            vars(1) = time - (1/DeltaT)*2*pow(time, 2) + (1/pow(DeltaT, 2))*pow(time, 3);           // x0dot coef
            return vars;
        } else if (poly_vars_.at(i-1).size() == 1) {
            // Initial time derivative is fixed at 0
            vector_t vars = vector_t::Zero(3);
            vars(0) = 1 - (1/pow(DeltaT, 2))*3*pow(time,2) + (1/pow(DeltaT, 3))*2*pow(time, 3);     // x0 coef
            vars(1) = (1/pow(DeltaT, 2))*3*pow(time,2) - (1/pow(DeltaT, 3))*2*pow(time, 3);         // x1 coef

            vars(2) = -(1/DeltaT)*pow(time, 2) + (1/pow(DeltaT, 2))*pow(time, 3);                   // x1dot coef

            return vars;
        } else {
            vector_t vars = vector_t::Zero(4);
            vars(0) = 1 - (1/pow(DeltaT, 2))*3*pow(time,2) + (1/pow(DeltaT, 3))*2*pow(time, 3);     // x0 coef
            vars(2) = (1/pow(DeltaT, 2))*3*pow(time,2) - (1/pow(DeltaT, 3))*2*pow(time, 3);         // x1 coef

            vars(1) = time - (1/DeltaT)*2*pow(time, 2) + (1/pow(DeltaT, 2))*pow(time, 3);           // x0dot coef
            vars(3) = -(1/DeltaT)*pow(time, 2) + (1/pow(DeltaT, 2))*pow(time, 3);                   // x1dot coef
            return vars;
        }

    }

    bool Spline::IsConstantPoly(int idx) const {
        return (poly_vars_.at(idx).size() == 1 && poly_vars_.at(idx-1).size() == 1);
    }

    std::pair<int, int> Spline::GetVarsIndexEnd(double time) const {
        // Given the index of the switching time, determine the start index of the (minimal) vars vector that
        // describes this polynomial

        if (time == poly_times_.at(0)) {
            int num_vars_affecting = 1;
            int vars_idx = 1;
            return std::pair<double, double>(vars_idx, num_vars_affecting);
        }

        int idx = GetPolyIdx(time);

        int vars_idx = 0;
        for (int i = 0; i <= idx; i++) {
            vars_idx += poly_vars_.at(i).size();
            if (i+1 < poly_vars_.size() && poly_vars_.at(i+1).size() == 1 && poly_vars_.at(i).size() == 1) {
                i++;
            }
        }

        int num_vars_affecting = 0;
        if (IsConstantPoly(idx)) {
            num_vars_affecting = 1;
        } else {
            num_vars_affecting = poly_vars_.at(idx).size() + poly_vars_.at(idx - 1).size();
        }

        return std::make_pair(vars_idx, num_vars_affecting);

    }

    int Spline::GetNumPolyVars(int idx) const {
        return poly_vars_.at(idx).size();
    }

    int Spline::GetNumPolyTimes() const {
        return poly_times_.size();
    }

    void Spline::UpdatePolyVar(int idx, const std::vector<double>& vars) {
        if (vars.size() != poly_vars_.at(idx).size()) {
            throw std::runtime_error("Number of variables does not match current number of variables.");
        }

        poly_vars_.at(idx) = vars;
    }

    const std::vector<std::vector<double>>& Spline::GetPolyVars() const {
        return poly_vars_;
    }

    const std::vector<double>& Spline::GetPolyTimes() const {
        return poly_times_;
    }

    int Spline::GetNumConstant() const {
        return num_constant_;
    }

    void Spline::SetAllPositions(double position) {
        for (int i = 0; i < poly_vars_.size(); i++) {
            poly_vars_.at(i).at(0) = position;
        }
    }

    int Spline::GetNumNonConstantValParams() const {
        int num = 0;
        for (const auto& poly : poly_vars_) {
            if (poly.size() != 1) {
                num += 1;
            }
        }

        return num;
    }

    double Spline::GetEndTime() const {
        return poly_times_.at(poly_times_.size()-1);
    }

    void Spline::AddPoly(double time) {
        int end_idx = poly_vars_.size()-1;
        std::vector<double> const ZERO_CONSTANT = {0};
        std::vector<double> const ZERO_POLY = {0, 0};

        if ((poly_vars_.at(end_idx).size() == 1 && poly_vars_.at(end_idx-1).size() != 1)) {
            poly_vars_.push_back(ZERO_CONSTANT);
            poly_times_.push_back(poly_times_.at(poly_times_.size()-1) + time);
        } else if (end_idx >= num_polys_ - 2 && poly_vars_.at(end_idx - num_polys_ + 2).size() == 2 &&
                    poly_vars_.at(end_idx).size() == 2) {
            poly_vars_.push_back(ZERO_CONSTANT);
            num_all_poly_vars++;
            num_constant_++;
            total_poly_++;
            poly_times_.push_back(poly_times_.at(poly_times_.size()-1) + time/num_polys_);
        } else {
            poly_vars_.push_back(ZERO_POLY);
            total_poly_++;
            num_all_poly_vars+=2;
            poly_times_.push_back(poly_times_.at(poly_times_.size()-1) + time/num_polys_);
        }
    }

    void Spline::RemoveUnused(double time) {
        for (int i = 0; i < poly_times_.size() - 1; i++) {
            if (poly_times_.at(i+1) < time && poly_times_.at(i) < time) {
                if (IsConstantPoly(i+1)) {
                    // Then we are getting rid of the first part of a constant, which is just an internal change

                } else if (poly_vars_.at(i).size() == 1) {
                    // Then we are removing the last part of a constant segment
                    total_poly_--;
                    num_constant_--;
                    num_all_poly_vars--;
                } else {
                    total_poly_--;
                    num_all_poly_vars -= 2;
                }
                poly_times_.erase(poly_times_.begin() + i);
                poly_vars_.erase(poly_vars_.begin() + i);
                i--;
            }
        }
    }

} // mpc