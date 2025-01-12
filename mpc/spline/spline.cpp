//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <cmath>
#include <cassert>
#include <iostream>

#include "spline/spline.h"

namespace mpc {
    Spline::Spline(int num_polys, const std::vector<double> &times, bool start_on_poly,
                   const SplineType& type) :
                    start_on_constant_(!start_on_poly), num_polys_(num_polys), type_(type){

        std::vector<double> const ZERO_CONSTANT = {0};
        std::vector<double> const ZERO_POLY = {0, 0};

        num_all_poly_vars = 0;
        total_poly_ = 0;
        num_constant_ = 0;

        poly_vars_.push_back(ZERO_CONSTANT);
        poly_times_.push_back(0);
        start_pair_ = start_on_constant_;
        end_pair_ = true;

        if (type == SplineType::Constants) {
            for (int idx = 0; idx < times.size(); idx++) {
                poly_vars_.push_back(ZERO_CONSTANT);
                num_constant_++;
                total_poly_++;
                poly_times_.push_back(times.at(idx));
                if (!(idx%2)){
                    num_all_poly_vars++;
                }
            }
            if (!(times.size()%2) && start_on_constant_){
                num_all_poly_vars++;
                end_pair_ = false;
            }
            num_constant_++;
            total_poly_++;
            if (!start_on_constant_) {
                num_all_poly_vars++;
            }
        } else {

            num_all_poly_vars++;
            num_constant_++;

            int mut_idx = 0;
            for (int time_idx = 0; time_idx < times.size(); time_idx++) {
                if (start_on_constant_) {
                    if ((time_idx) % 2 == 0) {      // On a constant
                        poly_vars_.push_back(ZERO_CONSTANT);
                        num_constant_++;
                        total_poly_++;
                        num_all_poly_vars++;

                        if (time_idx != 0) {
                            poly_times_.push_back(times.at(time_idx - 1));
                            poly_vars_.push_back(ZERO_CONSTANT);
                        }
                        poly_times_.push_back(times.at(time_idx));
                    } else {        // On a polynomial
                        total_poly_ += num_polys_;
                        for (int i = 0; i < num_polys_ - 1; i++) {
                            poly_vars_.push_back(ZERO_POLY);
                            poly_times_.push_back((times.at(time_idx) - times.at(time_idx - 1)) / num_polys
                                                  + poly_times_.at(poly_times_.size() - 1));

                            num_all_poly_vars += 2;
                        }

//                    poly_times_.push_back(times.at(time_idx));
                    }
                } else {
                    if ((time_idx + 1) % 2 == 0) {      // On a constant
                        poly_vars_.push_back(ZERO_CONSTANT);
                        poly_vars_.push_back(ZERO_CONSTANT);
                        num_constant_++;
                        total_poly_++;
                        total_poly_++;

                        num_all_poly_vars++;

                        poly_times_.push_back(times.at(time_idx - 1));
                        poly_times_.push_back(times.at(time_idx));
                    } else {        // On a polynomial
                        total_poly_ += num_polys_;

                        for (int i = 0; i < num_polys_ - 1; i++) {
                            if (time_idx == 0) {
                                poly_times_.push_back(
                                        times.at(time_idx) / num_polys + poly_times_.at(poly_times_.size() - 1));
                            } else {
                                poly_times_.push_back((times.at(time_idx) - times.at(time_idx - 1)) / num_polys
                                                      + poly_times_.at(poly_times_.size() - 1));
                            }

                            poly_vars_.push_back(ZERO_POLY);

                            num_all_poly_vars += 2;

                        }
                    }
                }
            }
            // Add final elements
            if ((start_on_constant_ && (times.size()) % 2 == 0) || (!start_on_constant_ && (times.size()-1) % 2 == 0)) {
                poly_vars_.push_back(ZERO_CONSTANT);
                poly_times_.push_back(times.at(times.size()-1));
                num_constant_++;
                num_all_poly_vars++;
            }
        }

        switch (type) {
            case PositionZ:
                throw std::runtime_error("Not supported yet");
                break;
            case Force:
                num_all_poly_vars -= num_constant_;
                for (const auto& poly_var : poly_vars_) {
                    if (poly_var.size() == 1) {
                        mut_flags_.emplace_back(false);
                    } else {
                        mut_flags_.emplace_back(true);
                    }
                }
                break;
            case Normal:
            case Constants:
            default:
                for (const auto& poly_var : poly_vars_) {
                    mut_flags_.emplace_back(true);
                }
                break;
        }
        assert(poly_times_.size() == poly_vars_.size() && poly_times_.size() == mut_flags_.size());
    }

    Spline& Spline::operator=(const Spline& spline) {
        poly_vars_ = spline.poly_vars_;
        poly_times_ = spline.poly_times_;
        mut_flags_ = spline.mut_flags_;
        total_poly_ = spline.total_poly_;
        start_on_constant_ = spline.start_on_constant_;
        num_constant_ = spline.num_constant_;
        num_polys_ = spline.num_polys_;
        num_all_poly_vars = spline.num_all_poly_vars;
        type_ = spline.type_;
        end_pair_ = spline.end_pair_;
        start_pair_ = spline.start_pair_;

        return *this;
    }

    double Spline::ValueAt(double time) const {
        std::array<double, POLY_ORDER> vars = GetPolyVars(time);

        int i = GetPolyIdx(time);
        double DeltaT = poly_times_.at(i) - poly_times_.at(i-1);
        time = time - poly_times_.at(i-1);

        return EvalPoly(vars, time, DeltaT);
    }

    void Spline::SetPolyVars(int idx, const std::vector<double>& vars) {
        assert(idx >= 0);
        assert(idx <= poly_times_.size());

        if (vars.size() != poly_vars_.at(idx).size()) {
            throw std::runtime_error("Provided polynomial does not match the order of the current polynomial.");
        }

        if (!IsMutable(idx)) {
            throw std::runtime_error("Chosen polynomial variables are not mutable.");
        }

        poly_vars_.at(idx) = vars;

        if (idx == 0) {
            // Don't need to check prev val
            if (poly_vars_.at(idx).size() == 1 && poly_vars_.at(idx + 1).size() == 1 && IsMutable(idx+1)) {
                poly_vars_.at(idx + 1) = vars;
            }
        } else if (idx == poly_times_.size()) {
            // Don't need to check after val
            if (poly_vars_.at(idx).size() == 1 && poly_vars_.at(idx - 1).size() == 1 && IsMutable(idx-1)) {
                poly_vars_.at(idx - 1) = vars;
            }
        } else {
            // Check for constants
            if (poly_vars_.at(idx).size() == 1 && poly_vars_.at(idx - 1).size() == 1 && IsMutable(idx-1)) {
                poly_vars_.at(idx - 1) = vars;
            } else if (idx + 1 < poly_vars_.size() && poly_vars_.at(idx).size() == 1 &&
                    poly_vars_.at(idx + 1).size() == 1 && IsMutable(idx)) {
                poly_vars_.at(idx + 1) = vars;
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

        std::cerr << "SetAllSplineVars not updated for mutability." << std::endl;

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
                if (mut_flags_.at(poly)) {
                    all_vars(idx) = poly_vars_.at(poly).at(j);
                    idx++;
                }
            }
            // Constants are always followed by redundant constants
            if (type_ == Constants) {
                if ((start_pair_ && poly == 0) || (end_pair_ && poly == poly_vars_.size()-2)) {
                    poly++;
                } else if (poly != 0 && poly != poly_vars_.size()-2) {
                    poly++;
                }
            } else {
                if (poly_vars_.at(poly).size() == 1 && poly < poly_vars_.size() - 1 &&
                    poly_vars_.at(poly + 1).size() == 1) {
                    poly++;
                }
            }
        }

        return all_vars;
    }

    int Spline::GetPolyIdx(double time) const {
        if (time - poly_times_.at(0) < 0 && time - poly_times_.at(0) > -1e-3) {
            time = poly_times_.at(0);
        } else if (time < poly_times_.at(0)) {
            throw std::runtime_error("Invalid time to query spline, time too small.");
        }

        if (time > poly_times_.at(poly_times_.size()-1) && (time - poly_times_.at(poly_times_.size()-1)) < 1e-3) {
            time = poly_times_.at(poly_times_.size()-1);
        }

        for (int i = 1; i < poly_times_.size(); i++) {     // Start at i = 1 because t = 0  at 0
            if (time < poly_times_.at(i)) {     // Note: removed equality here
                return i;
            }
        }

        if (time == poly_times_.back()) {
            return poly_times_.size() - 1;
        }

        // The spine is fully defined in the switching time, so if we got here then it's an invalid time.
        throw std::runtime_error("Invalid time to query the spline. Time too large.");
    }

    vector_t Spline::GetPolyVarsLin(double time) const {
        int i = GetPolyIdx(time);
        if (type_ == Constants) {
            int idx = i;
            if (start_pair_ && idx == 1) {
                vector_t vars(1);
                vars << 1;
                return vars;
            }

            if (end_pair_ && idx == poly_times_.size()-1) {
                vector_t vars(1);
                vars << 1;
                return vars;
            }

            if (start_pair_ && idx%2) {
                vector_t vars(1);
                vars << 1;
                return vars;
            }

            if (end_pair_ && !(idx%2)) {  // !start_pair
                vector_t vars(1);
                vars << 1;
                return vars;
            }

            if (start_pair_) {
                assert(!(idx%2));
                const double DeltaT = poly_times_.at(i) - poly_times_.at(i - 1);
                time = time - poly_times_.at(i - 1);

                vector_t vars(2);
                vars(0) = Getx0Coef(time, DeltaT);
                vars(1) = Getx1Coef(time, DeltaT);

                return vars;
            }

            if (end_pair_) {
                assert(idx%2);
                const double DeltaT = poly_times_.at(i) - poly_times_.at(i - 1);
                time = time - poly_times_.at(i - 1);

                vector_t vars(2);
                vars(0) = Getx0Coef(time, DeltaT);
                vars(1) = Getx1Coef(time, DeltaT);

                return vars;
            }

            if (idx%2) {
                const double DeltaT = poly_times_.at(i) - poly_times_.at(i - 1);
                time = time - poly_times_.at(i - 1);

                vector_t vars(2);
                vars(0) = Getx0Coef(time, DeltaT);
                vars(1) = Getx1Coef(time, DeltaT);

                return vars;
            } else {
                vector_t vars(1);
                vars << 1;
                return vars;
            }
        } else {
            if (IsConstantPoly(i)) {
                vector_t vars(1);
                vars << 1;
                return vars;
            }

            // Polynomial
            const double DeltaT = poly_times_.at(i) - poly_times_.at(i - 1);
            time = time - poly_times_.at(i - 1);

            if (mut_flags_.at(i) && mut_flags_.at(i - 1)) {
                if (poly_vars_.at(i).size() == 1) {
                    // Final time derivative is fixed at 0
                    vector_t vars = vector_t::Zero(3);
                    vars(0) = Getx0Coef(time, DeltaT);          // x0 coef
                    vars(2) = Getx1Coef(time, DeltaT);         // x1 coef
                    vars(1) = Getx0dotCoef(time, DeltaT);           // x0dot coef
                    return vars;
                } else if (poly_vars_.at(i - 1).size() == 1) {
                    // Initial time derivative is fixed at 0
                    vector_t vars = vector_t::Zero(3);
                    vars(0) = Getx0Coef(time, DeltaT);     // x0 coef
                    vars(1) = Getx1Coef(time, DeltaT);         // x1 coef
                    vars(2) = Getx1dotCoef(time, DeltaT);                   // x1dot coef

                    return vars;
                } else {
                    vector_t vars = vector_t::Zero(4);
                    vars(0) = Getx0Coef(time, DeltaT);     // x0 coef
                    vars(2) = Getx1Coef(time, DeltaT);         // x1 coef

                    vars(1) = Getx0dotCoef(time, DeltaT);           // x0dot coef
                    vars(3) = Getx1dotCoef(time, DeltaT);                   // x1dot coef
                    return vars;
                }
            } else if (mut_flags_.at(i)) {
                if (poly_vars_.at(i).size() == 1) {
                    // Final time derivative is fixed at 0
                    vector_t vars = vector_t::Zero(1);
                    vars(0) = Getx1Coef(time, DeltaT);         // x1 coef
                    return vars;
                } else if (poly_vars_.at(i - 1).size() == 1) {
                    // Initial time derivative is fixed at 0
                    vector_t vars = vector_t::Zero(2);
                    vars(0) = Getx1Coef(time, DeltaT);         // x1 coef
                    vars(1) = Getx1dotCoef(time, DeltaT);                   // x1dot coef
                    return vars;
                } else {
                    vector_t vars = vector_t::Zero(2);
                    vars(0) = Getx1Coef(time, DeltaT);         // x1 coef
                    vars(1) = Getx1dotCoef(time, DeltaT);                   // x1dot coef
                    return vars;
                }
            } else if (mut_flags_.at(i - 1)) {
                if (poly_vars_.at(i).size() == 1) {
                    // Final time derivative is fixed at 0
                    vector_t vars = vector_t::Zero(2);
                    vars(0) = Getx0Coef(time, DeltaT);     // x0 coef
                    vars(1) = Getx0dotCoef(time, DeltaT);           // x0dot coef
                    return vars;
                } else if (poly_vars_.at(i - 1).size() == 1) {
                    // Initial time derivative is fixed at 0
                    vector_t vars = vector_t::Zero(1);
                    vars(0) = Getx0Coef(time, DeltaT);     // x0 coef
                    return vars;
                } else {
                    vector_t vars = vector_t::Zero(2);
                    vars(0) = Getx0Coef(time, DeltaT);     // x0 coef
                    vars(1) = Getx0dotCoef(time, DeltaT);           // x0dot coef
                    return vars;
                }
            }
            throw std::runtime_error("Reached bad spot. PolyVarsLin");
        }
    }

    bool Spline::IsConstantPoly(int idx) const {
        return (poly_vars_.at(idx).size() == 1 && poly_vars_.at(idx-1).size() == 1);
        //&& poly_vars_.at(idx).at(0) == poly_vars_.at(idx-1).at(0));
    }

    std::pair<int, int> Spline::GetVarsIndexEnd(double time) const {
        // Given the index of the switching time, determine the start index of the (minimal) vars vector that
        // describes this polynomial

        int idx = GetPolyIdx(time);

        if (!mut_flags_.at(idx) && !mut_flags_.at(idx-1)) {
            throw std::runtime_error("The polynomial is not mutable at the provided time.");
        }

        if (type_ == Constants) {
            if (start_pair_ && idx == 1) {
                return std::make_pair(1,1);
            }

            if (end_pair_ && idx == poly_times_.size()-1) {
                return std::make_pair(std::floor(idx/2)+1, 1);
            }

            if (start_pair_ && idx%2) {
                return std::make_pair(std::floor(idx/2)+1, 1);
            }

            if (end_pair_ && !(idx%2)) { // !start_pair_
                return std::make_pair((idx-2)/2+2, 1);
            }

            if (start_pair_) {
                assert(!(idx%2));
                return std::make_pair(((idx-2)/2)+2, 2);
            }

            if (end_pair_) {
                assert(idx%2);
                return std::make_pair((idx-1)/2+2, 2);
            }

            if (idx%2) {
                return std::make_pair(2 + (idx-1)/2, 2);
            } else {
                return std::make_pair(2 + (idx-2)/2, 1);
            }
        }

        int vars_idx = 0;
        for (int i = 0; i <= idx; i++) {
            if (mut_flags_.at(i)) {
                vars_idx += poly_vars_.at(i).size();
                if (i + 1 < poly_vars_.size() && poly_vars_.at(i + 1).size() == 1 && poly_vars_.at(i).size() == 1) {
                    i++;
                }
            }
        }

        int num_vars_affecting = 0;
        if (IsConstantPoly(idx)) {
            num_vars_affecting = 1;
        } else {
            num_vars_affecting = static_cast<int>(mut_flags_.at(idx))*poly_vars_.at(idx).size() +
                                    static_cast<int>(mut_flags_.at(idx-1))*poly_vars_.at(idx - 1).size();
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
            if (mut_flags_.at(i)) {
                poly_vars_.at(i).at(0) = position;
            }
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
        bool added_constant = false;

        if (type_ == Constants) {
            poly_vars_.emplace_back(poly_vars_.at(poly_vars_.size()-1));
            total_poly_++;
            num_constant_++;
            if (start_pair_ && !(poly_times_.size()%2) || !start_pair_ && (poly_times_.size()%2)) {
                num_all_poly_vars++;
                end_pair_ = false;
            } else {
                end_pair_ = true;
            }
            poly_times_.push_back(poly_times_.at(poly_times_.size()-1) + time);
            mut_flags_.push_back(true);
        } else {
            if (poly_vars_.at(end_idx).size() == 1 && poly_vars_.at(end_idx - 1).size() == 1) {
                // External and internal change
                for (int i = 0; i < num_polys_ - 1; i++) {
                    poly_vars_.push_back(ZERO_POLY);
                    total_poly_++;
                    num_all_poly_vars += 2;
                    poly_times_.push_back(poly_times_.at(poly_times_.size() - 1) + time / num_polys_);
                    mut_flags_.emplace_back(true);
                }

                // Internal and external change
                poly_vars_.push_back(ZERO_CONSTANT);
                poly_times_.push_back(poly_times_.at(poly_times_.size() - 1) + time / num_polys_);
                num_constant_++;
                total_poly_++;

                if (type_ == Normal) {
                    mut_flags_.emplace_back(true);
                    num_all_poly_vars++;
                } else {
                    mut_flags_.emplace_back(false);
                }
            } else {
                // Internal and external change
                poly_vars_.push_back(ZERO_CONSTANT);
                poly_times_.push_back(poly_times_.at(poly_times_.size() - 1) + time);

                if (type_ == Normal) {
                    mut_flags_.emplace_back(true);
                } else {
                    mut_flags_.emplace_back(false);
                }
            }
        }
    }

    void Spline::RemoveUnused(double time) {
        for (int i = 0; i < poly_times_.size() - 1; i++) {
            if (poly_times_.at(i+1) <= time && poly_times_.at(i) < time) { // make the second inequality also equality
                if (IsConstantPoly(i+1)) {
                    // Then we are getting rid of the first part of a constant, which is just an internal change
                    if (!start_pair_ && type_ == Constants) {
                        num_all_poly_vars--;
                        num_constant_--;
                        total_poly_--;
                        start_pair_ = true;
                    } else {
                        start_pair_ = false;
                    }
                } else if (poly_vars_.at(i).size() == 1) {
                    // Then we are removing the last part of a constant segment
                    total_poly_--;
                    num_constant_--;
                    if (type_ == Normal || type_ == Constants) {
                        num_all_poly_vars--;
                    }
                } else {
                    total_poly_--;
                    num_all_poly_vars -= 2;
                }
                poly_times_.erase(poly_times_.begin() + i);
                poly_vars_.erase(poly_vars_.begin() + i);
                mut_flags_.erase(mut_flags_.begin() + i);
                i--;
            }
        }
    }

    bool Spline::IsConstant(double time) const {
        int idx = GetPolyIdx(time);
        return IsConstantPoly(idx) || (time == poly_times_.at(0) && poly_vars_.at(0).size() == 1);
    }

    bool Spline::IsMutable(int idx) const {
        return mut_flags_.at(idx);
    }

    double Spline::ComputePartialWrtTime(double time, int time_idx) const {
        const int idx = GetPolyIdx(time);
        const double time_spline = time - poly_times_.at(idx-1);
        const double DeltaT = poly_times_.at(idx) - poly_times_.at(idx-1);

        if (time_idx == idx) {
            const double x0 = poly_vars_.at(idx-1).at(0);
            const double x1 = poly_vars_.at(idx).at(0);
            double x0dot = 0;
            double x1dot = 0;
            if (poly_vars_.at(idx-1).size() == 2) {
                x0dot = poly_vars_.at(idx-1).at(1);
            }
            if (poly_vars_.at(idx).size() == 2) {
                x1dot = poly_vars_.at(idx).at(1);
            }

            const double da2dt2 = 2*pow(DeltaT, -3)*(3*(x0 - x1) + DeltaT*(2*x0dot - x1dot)) -
                                  pow(DeltaT, -2)*(3*(x0 - x1) + (2*x0dot + x1dot));
            const double da3dt2 = -3*pow(DeltaT, -4)*(2*(x0 - x1) + DeltaT*(x0dot + x1dot)) +
                                  pow(DeltaT, -3)*(2*(x0 - x1) + (x0dot + x1dot));

            return da2dt2*pow(time_spline,2) + da3dt2* pow(time_spline, 3);

        } else if (time_idx == idx - 1) {
            const double x0 = poly_vars_.at(idx-1).at(0);
            const double x1 = poly_vars_.at(idx).at(0);
            double x0dot = 0;
            double x1dot = 0;
            if (poly_vars_.at(idx-1).size() == 2) {
                x0dot = poly_vars_.at(idx-1).at(1);
            }
            if (poly_vars_.at(idx).size() == 2) {
                x1dot = poly_vars_.at(idx).at(1);
            }

            const double da2dt1 = -2*pow(DeltaT, -3)*(3*(x0 - x1) + DeltaT*(2*x0dot - x1dot)) -
                                  pow(DeltaT, -2)*(3*(x0 - x1) - (2*x0dot + x1dot));
            const double da3dt1 = 3*pow(DeltaT, -4)*(2*(x0 - x1) + DeltaT*(x0dot + x1dot)) +
                                  pow(DeltaT, -3)*(2*(x0 - x1) - (x0dot + x1dot));

            return da2dt1*pow(time_spline,2) + da3dt1* pow(time_spline, 3);

        } else {
            return 0;
        }
    }

    vector_t Spline::ComputeCoefPartialWrtTime(double time, int time_idx) const {
        const int idx = GetPolyIdx(time);

        bool wrt_end = false;
        if (idx == time_idx) {
            wrt_end = true;
        }

        const double time_spline = time - poly_times_.at(idx-1);
        const double DeltaT = poly_times_.at(idx) - poly_times_.at(idx-1);

        int vars_index, vars_affecting;
        std::tie(vars_index, vars_affecting) = GetVarsIndexEnd(time);

        vector_t coef_partials(vars_affecting);
        if (time_idx != idx || time_idx != idx-1) {
            coef_partials.setZero();
            return coef_partials;
        }

        coef_partials(0) = Getx0CoefPartial(time, DeltaT, !wrt_end);
        // Now need to determine if it is x0dot, or x1 next
        if (poly_vars_.at(idx-1).size() == 2) {
            coef_partials(1) = Getx0dotCoefPartial(time, DeltaT, !wrt_end);
            coef_partials(2) = Getx1CoefPartial(time, DeltaT, wrt_end);
            if (poly_vars_.at(idx).size() == 2) {
                coef_partials(3) = Getx1dotCoefPartial(time, DeltaT, wrt_end);
            }
        } else {
            coef_partials(1) = Getx1CoefPartial(time, DeltaT, wrt_end);
            if (poly_vars_.at(idx).size() == 2) {
                coef_partials(2) = Getx1dotCoefPartial(time, DeltaT, wrt_end);
            }
        }

        return coef_partials;
    }

    double Spline::Getx0Coef(double time, double DeltaT) const {
        return 1 - (1 / pow(DeltaT, 2)) * 3 * pow(time, 2) +
        (1 / pow(DeltaT, 3)) * 2 * pow(time, 3);
    }

    double Spline::Getx1Coef(double time, double DeltaT) const {
        return (1 / pow(DeltaT, 2)) * 3 * pow(time, 2) -
        (1 / pow(DeltaT, 3)) * 2 * pow(time, 3);
    }

    double Spline::Getx0dotCoef(double time, double DeltaT) const {
        return time - (1 / DeltaT) * 2 * pow(time, 2) +
        (1 / pow(DeltaT, 2)) * pow(time, 3);
    }

    double Spline::Getx1dotCoef(double time, double DeltaT) const {
        return -(1 / DeltaT) * pow(time, 2) +
        (1 / pow(DeltaT, 2)) * pow(time, 3);
    }

    double Spline::Getx0CoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return -6*pow(DeltaT, -3)* pow(time, 2) + 6*pow(DeltaT, -2)*time
            - 6*pow(DeltaT, 2)*pow(time,3) - 6*pow(DeltaT, 3)* pow(time, 2);
        } else {
            return 6*pow(DeltaT, -3)* pow(time, 2) + 6* pow(DeltaT, 2)* pow(time, 3);
        }
    }

    double Spline::Getx1CoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return 6*pow(DeltaT, -3)* pow(time, 2) - 6* pow(DeltaT, -2)*time
                    - 6* pow(DeltaT, -4)* pow(time, 3) + 6*pow(DeltaT, -3)* pow(time, 2);
        } else {
            return -6* pow(DeltaT, -3)* pow(time, 2) + 6* pow(DeltaT, -4)* pow(time, 3);
        }
    }

    double Spline::Getx0dotCoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return -1 -2*pow(DeltaT,-2)*pow(time, 2) + 4* pow(DeltaT, -1)*time
                + 2* pow(DeltaT, -3)* pow(time, 3) - 3*pow(DeltaT, -2)*pow(time, 2);
        } else {
            return 2*pow(DeltaT, -2)* pow(time, 2) - 2* pow(DeltaT, -3)* pow(time, 3);
        }
    }

    double Spline::Getx1dotCoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return -pow(DeltaT, -2)* pow(time, 2) + 2*pow(DeltaT, -1)*time
                + 2*pow(DeltaT,-3)*pow(time, 3) - 3*pow(DeltaT, -2)*pow(time, 2);
        } else {
            return -pow(DeltaT, -2)* pow(time, 2) - 2*pow(DeltaT, -3)*pow(time, 3);
        }
    }

    bool Spline::IsStartPairConstant() const {
        return start_pair_;
    }

    bool Spline::IsEndPairConstant() const {
        return end_pair_;
    }

    void Spline::ChangeSplineTimes(const std::vector<double>& times) {
        // Need to assign the times at every point associated with a constant
        int times_idx = 0;
        for (int i = 0; i < poly_times_.size(); i++) {
            if (poly_vars_.at(i).size() == 1) {
                poly_times_.at(i) = times.at(times_idx);
                times_idx++;
            } else if (times_idx > 0) {
                poly_times_.at(i) = poly_times_.at(i-1) + (times.at(times_idx) - times.at(times_idx-1))/num_polys_;
            } else {
                poly_times_.at(i + 1) = poly_times_.at(i) + (times.at(times_idx) - poly_times_.at(0))/num_polys_;
                i++;
            }
        }
    }

} // mpc