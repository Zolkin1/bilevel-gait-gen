//
// Created by zolkin on 1/30/24.
//
#include <assert.h>

#include <iostream>

#include "end_effector_splines.h"

namespace mpc {

    SplineTimes::SplineTimes(double time, TimeType type) {
        time_ = time;
        type_ = type;
    }

    SplineTimes::SplineTimes() {
        time_ = -1;
        type_ = Inter;
    }

    double SplineTimes::GetTime() const {
        return time_;
    }

    TimeType SplineTimes::GetType() const {
        return type_;
    }

    void SplineTimes::SetTime(double time) {
        time_ = time;
    }

    EndEffectorSplines::EndEffectorSplines(int num_contacts,
                                           const std::vector<double>& times,
                                           bool start_in_contact,
                                           int num_force_polys)
                                           : num_force_polys_(num_force_polys) {
        const vector_2t ZERO_VARS = vector_2t::Zero();

        if (num_force_polys_ < 2) {
            throw std::runtime_error("The number of force polynomials between constant sections must be at least 2.");
        }

        int num_times_in_pattern = num_force_polys + 1;
        spline_stride_ = num_force_polys;
        if (num_force_polys % 2) {
            num_times_in_pattern++;
            spline_stride_ = num_force_polys;
        }

        assert(num_times_in_pattern % 2);

        for (int i = 0; i < num_times_in_pattern; i++) {
            if (!start_in_contact) {
                if (i == 0) {
                    force_type_pattern_.push_back(NoDeriv);
                    position_type_pattern_.push_back(NoDeriv);
                    z_position_type_pattern_.push_back(NoDeriv);
                    time_type_pattern_.push_back(LiftOff);
                } else if (i == 1) {
                    force_type_pattern_.push_back(Empty);
                    position_type_pattern_.push_back(Empty);
                    z_position_type_pattern_.push_back(FullDeriv);
                    time_type_pattern_.push_back(Inter);
                } else if (i == 2) {
                    force_type_pattern_.push_back(NoDeriv);
                    position_type_pattern_.push_back(NoDeriv);
                    z_position_type_pattern_.push_back(NoDeriv);
                    time_type_pattern_.push_back(TouchDown);
                } else {
                    force_type_pattern_.push_back(FullDeriv);
                    position_type_pattern_.push_back(Empty);
                    z_position_type_pattern_.push_back(Empty);
                    time_type_pattern_.push_back(Inter);
                }
            } else {
                if (i == 0) {
                    force_type_pattern_.push_back(NoDeriv);
                    position_type_pattern_.push_back(NoDeriv);
                    z_position_type_pattern_.push_back(NoDeriv);
                    time_type_pattern_.push_back(TouchDown);
                } else if (i < num_force_polys) {
                    force_type_pattern_.push_back(FullDeriv);
                    position_type_pattern_.push_back(Empty);
                    z_position_type_pattern_.push_back(Empty);
                    time_type_pattern_.push_back(Inter);
                } else if ( i == num_times_in_pattern - 1) {
                    force_type_pattern_.push_back(Empty);
                    position_type_pattern_.push_back(Empty);
                    z_position_type_pattern_.push_back(FullDeriv);
                    time_type_pattern_.push_back(Inter);
                } else {
                    force_type_pattern_.push_back(NoDeriv);
                    position_type_pattern_.push_back(NoDeriv);
                    z_position_type_pattern_.push_back(NoDeriv);
                    time_type_pattern_.push_back(LiftOff);
                }
            }
        }


        for (int coord = 0; coord < POS_VARS; coord++) {
            int i = 0;
            int j = 0;
            int k = 1;
            if (coord < 2) {
                while (i < num_contacts) {
                    forces_.at(coord).emplace_back(force_type_pattern_.at(j % force_type_pattern_.size()), j, ZERO_VARS);
                    positions_.at(coord).emplace_back(position_type_pattern_.at(j % position_type_pattern_.size()), j,
                                                      ZERO_VARS);

                    if (force_type_pattern_.at(j % position_type_pattern_.size()) == FullDeriv) {
                        if (coord == 0) {
                            times_.emplace_back(
                                    times.at(i-1) + k * (times.at(i) - times.at(i - 1)) / (num_force_polys),
                                    time_type_pattern_.at(j % time_type_pattern_.size()));
                            k++;
                        }
                    } else if (force_type_pattern_.at(j % position_type_pattern_.size()) == Empty) {
                        if (coord == 0) {
                            times_.emplace_back(times.at(i - 1) + (times.at(i) - times.at(i - 1)) / 2,
                                             time_type_pattern_.at(j % time_type_pattern_.size()));
                        }
                    } else {
                        if (coord == 0) {
                            times_.emplace_back(times.at(i),
                                             time_type_pattern_.at(j % time_type_pattern_.size()));
                        }
                        i++;
                        k = 1;
                    }
                    j++;
                }
            } else {
                while (i < num_contacts) {
                    forces_.at(coord).emplace_back(force_type_pattern_.at(j % force_type_pattern_.size()), j, ZERO_VARS);
                    positions_.at(coord).emplace_back(z_position_type_pattern_.at(j % z_position_type_pattern_.size()), j,
                                                      ZERO_VARS);

                    if (force_type_pattern_.at(j % position_type_pattern_.size()) == FullDeriv ||
                        force_type_pattern_.at(j % position_type_pattern_.size()) == Empty) {
                    } else {
                        i++;
                    }
                    j++;
                }
            }
        }

        assert(forces_.at(1).size() == positions_.at(1).size());
        assert(forces_.at(1).size() == times_.size());
    }

    EndEffectorSplines& EndEffectorSplines::operator=(const EndEffectorSplines& ee_spline) {
        if (this == &ee_spline) {
            return *this;
        }

        this->times_ = ee_spline.times_;
        this->forces_ = ee_spline.forces_;
        this->positions_ = ee_spline.positions_;
        this->num_force_polys_ = ee_spline.num_force_polys_;
        this->spline_stride_ = ee_spline.spline_stride_;

        return *this;
    }

    double EndEffectorSplines::ValueAt(SplineType type, int coord, double time) const {
        const node_v& spline = SelectSpline(type, coord);
        const int lower_node = GetLowerNodeIdx(type, coord, time);
        const int upper_node = GetUpperNodeIdx(type, coord, time);

        if (upper_node == lower_node) {
            return spline.at(lower_node).GetVars()(0);
        }

        const double deltat = times_.at(upper_node).GetTime() - times_.at(lower_node).GetTime();
        const double time_spline = time - times_.at(lower_node).GetTime();
        const vector_4t vars {spline.at(lower_node).GetVars()(0),         // x0
                              spline.at(upper_node).GetVars()(0),         // x1
                              spline.at(lower_node).GetVars()(1),         // x0dot
                              spline.at(upper_node).GetVars()(1)};       // x1dot

        const double a2 = -(1 / pow(deltat, 2)) * 3 * (vars(0) - vars(1)) -
                    (1 / deltat) * (2 * vars(2) + vars(3));
        const double a3 = (1 / pow(deltat, 3)) * 2 * (vars(0) - vars(1)) +
                    (1 / pow(deltat, 2)) * (vars(2) + vars(3));

        const double val = vars(0) + vars(2) * time_spline + a2 * pow(time_spline, 2) + a3 * pow(time_spline, 3);

        return val;

    }

    vector_t EndEffectorSplines::GetPolyVarsLin(SplineType type, int coord, double time) const {
        const node_v& spline = SelectSpline(type, coord);
        const int lower_node = GetLowerNodeIdx(type, coord, time);
        const int upper_node = GetUpperNodeIdx(type, coord, time);

        if (lower_node == upper_node) {
            vector_t vars_lin(1);
            vars_lin(0) = 1;
            return vars_lin;
        }

        const double poly_time = time - times_.at(lower_node).GetTime();
        const double deltat = times_.at(upper_node).GetTime() - times_.at(lower_node).GetTime();

        vector_t vars_lin;
        if (type == Force) {
            if (spline.at(lower_node).GetType() == NoDeriv && spline.at(upper_node).GetType() == NoDeriv) {
                throw std::runtime_error("There is no mutable variables at the provided time.");
            }

            if (spline.at(lower_node).GetType() == NoDeriv && spline.at(upper_node).GetType() == FullDeriv) {
                vars_lin.resize(2);
                vars_lin(0) = Getx1Coef(poly_time, deltat);
                vars_lin(1) = Getx1dotCoef(poly_time, deltat);

                return vars_lin;

            } else if (spline.at(lower_node).GetType() == FullDeriv && spline.at(upper_node).GetType() == NoDeriv) {

                vars_lin.resize(2);
                vars_lin(0) = Getx0Coef(poly_time, deltat);
                vars_lin(1) = Getx0dotCoef(poly_time, deltat);

                return vars_lin;
            }

            vars_lin.resize(4);
            vars_lin(0) = Getx0Coef(poly_time, deltat);
            vars_lin(1) = Getx0dotCoef(poly_time, deltat);
            vars_lin(2) = Getx1Coef(poly_time, deltat);
            vars_lin(3) = Getx1dotCoef(poly_time, deltat);

            return vars_lin;
        } else {
            if (coord != 2) {
                if (forces_.at(coord).at(lower_node).GetType() == NoDeriv &&
                    forces_.at(coord).at(lower_node + 2).GetType() == NoDeriv) {
                    vars_lin.resize(2);
                    vars_lin(0) = Getx0Coef(poly_time, deltat);
                    vars_lin(1) = Getx1Coef(poly_time, deltat);

                    return vars_lin;
                } else {
                    vars_lin.resize(1);
                    vars_lin << 1;
                    return vars_lin;
                }
            } else {
                if (positions_.at(coord).at(lower_node).GetType() == NoDeriv &&
                    positions_.at(coord).at(upper_node).GetType() == FullDeriv) {
                    vars_lin.resize(3);
                    vars_lin(0) = Getx0Coef(poly_time, deltat);
                    vars_lin(1) = Getx1Coef(poly_time, deltat);
                    vars_lin(2) = Getx1dotCoef(poly_time, deltat);

                    return vars_lin;
                } else if (positions_.at(coord).at(lower_node).GetType() == FullDeriv &&
                            positions_.at(coord).at(upper_node).GetType() == NoDeriv) {
                    vars_lin.resize(3);
                    vars_lin(0) = Getx0Coef(poly_time, deltat);
                    vars_lin(1) = Getx0dotCoef(poly_time, deltat);
                    vars_lin(2) = Getx1Coef(poly_time, deltat);

                    return vars_lin;
                } else {
                    vars_lin.resize(1);
                    vars_lin << 1;
                    return vars_lin;
                }
            }
        }
    }

    std::pair<int, int> EndEffectorSplines::GetVarsIdx(SplineType type, int coord, double time) const {
        const node_v& spline = SelectSpline(type, coord);

        const int lower_node = GetLowerNodeIdx(type, coord, time);
        const int upper_node = GetUpperNodeIdx(type, coord, time);

        const std::vector<int> mut_nodes = GetMutableNodes(type, coord);

        int vars_idx = 0;
        if (type == Force) {
            for (int i = 0; i < mut_nodes.size(); i++) {
                if (mut_nodes.at(i) < lower_node) {
                    vars_idx = 2*(i+1);
                }
            }

            if (spline.at(lower_node).GetType() == NoDeriv && spline.at(upper_node).GetType() == NoDeriv) {
                throw std::runtime_error("There is no mutable variables at the provided time.");
            }

            if ((spline.at(lower_node).GetType() == NoDeriv && spline.at(upper_node).GetType() == FullDeriv) ||
                (spline.at(lower_node).GetType() == FullDeriv && spline.at(upper_node).GetType() == NoDeriv)) {
                return std::make_pair(vars_idx, 2);
            }

            if (lower_node == upper_node) {
                return std::make_pair(vars_idx, 1);
            }

            return std::make_pair(vars_idx, 4);
        } else {
            vars_idx--;
            for (int i = 0; i < mut_nodes.size(); i++) {
                if (mut_nodes.at(i) <= lower_node) {
                    vars_idx++;
                }
            }

            if (lower_node == upper_node) {
                return std::make_pair(vars_idx, 1);
            }

            if (coord != 2) {
                if (forces_.at(coord).at(lower_node).GetType() == NoDeriv &&
                    forces_.at(coord).at(lower_node + 2).GetType() == NoDeriv) {
                    return std::make_pair(vars_idx, 2);
                } else {
                    return std::make_pair(vars_idx, 1);
                }
            } else {
                for (int i = 0; i < mut_nodes.size(); i++) {
                    if (positions_.at(coord).at(mut_nodes.at(i)).GetType() == FullDeriv && mut_nodes.at(i) < lower_node) {
                        vars_idx++;
                    }
                }

                if ((positions_.at(coord).at(lower_node).GetType() == NoDeriv &&
                    positions_.at(coord).at(upper_node).GetType() == FullDeriv) ||
                        (positions_.at(coord).at(lower_node).GetType() == FullDeriv &&
                         positions_.at(coord).at(upper_node).GetType() == NoDeriv)) {
                    return std::make_pair(vars_idx, 3);
                } else if (forces_.at(coord).at(lower_node).GetType() == NoDeriv &&
                           forces_.at(coord).at(lower_node + 2).GetType() == NoDeriv) {
                    return std::make_pair(vars_idx, 2);
                } else {
                    return std::make_pair(vars_idx, 1);
                }
            }
        }

    }

    bool EndEffectorSplines::IsForceMutable(double time) const {
        const int lower_node = GetLowerNodeIdx(Force, 0, time);
        const int upper_node = GetUpperNodeIdx(Force, 0, time);
        if ((forces_.at(0).at(lower_node).GetType() == NoDeriv && forces_.at(0).at(upper_node).GetType() == NoDeriv)) {
            return false;
        } else {
            return true;
        }
    }

    void EndEffectorSplines::AddPoly(double additional_time) {
//        std::cerr << "poly added" << std::endl;

        const int num_nodes = GetNumNodes();
        const vector_2t ZERO_VEC = {0, 0};

        if (forces_.at(0).at(num_nodes-1).GetType() == NoDeriv && forces_.at(0).at(num_nodes-2).GetType() == FullDeriv) {
            for (int i = 0; i < 2; i++) {
                for (int coord = 0; coord < POS_VARS; coord++) {
                    if (i == 0) {
                        forces_.at(coord).emplace_back(Empty,
                                                       forces_.at(coord).at(forces_.at(coord).size() - 1).GetNodeIdx() +
                                                       1,
                                                       ZERO_VEC);

                        if (coord == 0) {
                            times_.emplace_back(times_.at(times_.size() - 1).GetTime() + additional_time / 2,
                                                Inter);
                        }

                        if (coord == 2) {
                            positions_.at(coord).emplace_back(FullDeriv,
                                                              positions_.at(coord).at(
                                                                      positions_.at(coord).size() - 1).GetNodeIdx() + 1,
                                                              ZERO_VEC);
                        } else {
                            positions_.at(coord).emplace_back(Empty,
                                                              positions_.at(coord).at(
                                                                      positions_.at(coord).size() - 1).GetNodeIdx() + 1,
                                                              ZERO_VEC);
                        }

                    } else {
                        forces_.at(coord).emplace_back(NoDeriv,
                                                       forces_.at(coord).at(forces_.at(coord).size() - 1).GetNodeIdx() +
                                                       1,
                                                       ZERO_VEC);
                        positions_.at(coord).emplace_back(NoDeriv,
                                                          positions_.at(coord).at(
                                                                  positions_.at(coord).size() - 1).GetNodeIdx() + 1,
                                                          ZERO_VEC);

                        if (coord == 0) {
                            times_.emplace_back(times_.at(times_.size() - 1).GetTime() + additional_time / 2,
                                                TouchDown);
                        }
                    }
                }
            }

        } else {
            for (int coord = 0; coord < POS_VARS; coord++) {
                for (int i = 0; i < num_force_polys_-1; i++) {
                    forces_.at(coord).emplace_back(FullDeriv,
                                                   forces_.at(coord).at(forces_.at(coord).size() - 1).GetNodeIdx() +
                                                   1,
                                                   ZERO_VEC);
                    positions_.at(coord).emplace_back(Empty,
                                                      positions_.at(coord).at(positions_.at(coord).size() - 1).GetNodeIdx() + 1,
                                                          ZERO_VEC);

                    if (coord == 0) {
                        times_.emplace_back(times_.at(times_.size()-1).GetTime() + additional_time/num_force_polys_,
                                            Inter);
                    }
                }
                forces_.at(coord).emplace_back(NoDeriv,
                                               forces_.at(coord).at(forces_.at(coord).size() - 1).GetNodeIdx() +
                                               1,
                                               ZERO_VEC);
                positions_.at(coord).emplace_back(NoDeriv,
                                                  positions_.at(coord).at(positions_.at(coord).size() - 1).GetNodeIdx() + 1,
                                                  ZERO_VEC);
                if (coord == 0) {
                    times_.emplace_back(times_.at(times_.size() - 1).GetTime() + additional_time / num_force_polys_,
                                     LiftOff);
                }
            }
        }

        for (int i = 1; i < times_.size(); i++) {
            assert(times_.at(i-1).GetTime() <= times_.at(i).GetTime());
        }
    }

    void EndEffectorSplines::RemovePoly(double start_time) {
        for (int i = 0; i < times_.size(); i++) {
            std::vector<int> mut_nodes = GetMutableNodes(Position, 0);

            for (int node = 1; node < mut_nodes.size(); node++) {
                if (times_.at(mut_nodes.at(node)).GetTime() < start_time) {
                    for (int j = 0; j < mut_nodes.at(node); j++) {
                        times_.erase(times_.begin());
                        for (int coord = 0; coord < POS_VARS; coord++) {
                            forces_.at(coord).erase(forces_.at(coord).begin());
                            positions_.at(coord).erase(positions_.at(coord).begin());
                        }
                    }
//                    std::cerr << "nodes removed" << std::endl;
                }
                mut_nodes = GetMutableNodes(Position, 0);
            }
        }

        for (int coord = 0; coord < POS_VARS; coord++) {
            assert(forces_.at(coord).at(0).GetType() != Empty);
            assert(positions_.at(coord).at(0).GetType() != Empty);
        }

//        for (int node = 1; node < mut_nodes.size(); node++) {
//            if (times_.at(mut_nodes.at(node-1)) < start_time && times_.at(mut_nodes.at(node)) < start_time) {
//                times_.erase(times_.begin() + node-1);
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    forces_.at(coord).erase(forces_.at(coord).begin() + node - 1);
//                    positions_.at(coord).erase(positions_.at(coord).begin() + node - 1);
//                }
//                node--;
//            }
//        }

//        for (int i = 1; i < times_.size(); i++) {
//            if (times_.at(i-1) < start_time && times_.at(i) < start_time) {
//                times_.erase(times_.begin() + i-1);
//                for (int coord = 0; coord < POS_VARS; coord++) {
//                    forces_.at(coord).erase(forces_.at(coord).begin() + i - 1);
//                    positions_.at(coord).erase(positions_.at(coord).begin() + i - 1);
//                }
//                i--;
//            }
//        }
    }

    double EndEffectorSplines::ComputePartialWrtTime(SplineType type, int coord, double time, int time_idx) const {
        const node_v& spline = SelectSpline(type, coord);
        const int upper_node = GetUpperNodeIdx(type, coord, time);
        const int lower_node = GetLowerNodeIdx(type, coord, time);

        const double deltat = times_.at(upper_node).GetTime() - times_.at(lower_node).GetTime();
        assert(deltat != 0);
        const double time_spline = time - times_.at(lower_node).GetTime();

        // Determine if this is related to the upper or lower index (or neither)
        // If this is a position then it is always a direct dependence
        //      (i.e. the time_idx corresponds with a node) or not at all
        // If this is a force then it can either be a direct dependence or an indirect
        //      (the inernal nodes shift based on the location of the contact points), or nothing

        // Start by determining if it is a direct dependence
        bool direct_dep = false;
        bool wrt_lower = false;

        if (type == Position) {
            direct_dep = true;

            const int node = ConvertContactNodeToSplineNode(time_idx);
            if (node == lower_node) {
                wrt_lower = true;
            }
        } else {
            const int node = ConvertContactNodeToSplineNode(time_idx);
            if (forces_.at(coord).at(node).GetType() == NoDeriv) {
                direct_dep = true;
            }

            if (direct_dep && node == lower_node) {
                wrt_lower = true;
            }
        }

        if (direct_dep && wrt_lower) {
            const double x0 = spline.at(lower_node).GetVars()(0);
            const double x1 = spline.at(upper_node).GetVars()(0);
            double x0dot = 0;
            double x1dot = 0;
            if (spline.at(lower_node).GetType() == FullDeriv) {
                x0dot = spline.at(lower_node).GetVars()(1);
            }
            if (spline.at(upper_node).GetType() == FullDeriv) {
                x1dot = spline.at(upper_node).GetVars()(1);
            }

            const double da2dt1 = -2*pow(deltat, -3)*(3*(x0 - x1) + deltat*(2*x0dot - x1dot)) -
                                  pow(deltat, -2)*(3*(x0 - x1) - (2*x0dot + x1dot));
            const double da3dt1 = 3*pow(deltat, -4)*(2*(x0 - x1) + deltat*(x0dot + x1dot)) +
                                  pow(deltat, -3)*(2*(x0 - x1) - (x0dot + x1dot));

            return da2dt1*pow(time_spline,2) + da3dt1* pow(time_spline, 3);
        }

        if (direct_dep && !wrt_lower) {
            const double x0 = spline.at(lower_node).GetVars()(0);
            const double x1 = spline.at(upper_node).GetVars()(0);
            double x0dot = 0;
            double x1dot = 0;
            if (spline.at(lower_node).GetType() == FullDeriv) {
                x0dot = spline.at(lower_node).GetVars()(1);
            }
            if (spline.at(upper_node).GetType() == FullDeriv) {
                x1dot = spline.at(upper_node).GetVars()(1);
            }

            const double da2dt2 = 2*pow(deltat, -3)*(3*(x0 - x1) + deltat*(2*x0dot - x1dot)) -
                                  pow(deltat, -2)*(3*(x0 - x1) + (2*x0dot + x1dot));
            const double da3dt2 = -3*pow(deltat, -4)*(2*(x0 - x1) + deltat*(x0dot + x1dot)) +
                                  pow(deltat, -3)*(2*(x0 - x1) + (x0dot + x1dot));

            return da2dt2*pow(time_spline,2) + da3dt2* pow(time_spline, 3);
        }

        // Not a direct dependency
        // TODO: Check this math!
        const int node = ConvertContactNodeToSplineNode(time_idx);
        if (forces_.at(coord).at(upper_node).GetType() != NoDeriv && node == upper_node + 1) {
            const double x0 = spline.at(lower_node).GetVars()(0);
            const double x1 = spline.at(upper_node).GetVars()(0);
            double x0dot = 0;
            double x1dot = 0;
            if (spline.at(lower_node).GetType() == FullDeriv) {
                x0dot = spline.at(lower_node).GetVars()(1);
            }
            if (spline.at(upper_node).GetType() == FullDeriv) {
                x1dot = spline.at(upper_node).GetVars()(1);
            }

            const double da2dt2 = 2*pow(deltat, -3)*(3*(x0 - x1) + deltat*(2*x0dot - x1dot)) -
                                  pow(deltat, -2)*(3*(x0 - x1) + (2*x0dot + x1dot));
            const double da3dt2 = -3*pow(deltat, -4)*(2*(x0 - x1) + deltat*(x0dot + x1dot)) +
                                  pow(deltat, -3)*(2*(x0 - x1) + (x0dot + x1dot));

            return da2dt2*pow(time_spline,2) + da3dt2* pow(time_spline, 3)/static_cast<double>(num_force_polys_);
        }

        if (forces_.at(coord).at(lower_node).GetType() != NoDeriv && node == lower_node - 1) {
            const double x0 = spline.at(lower_node).GetVars()(0);
            const double x1 = spline.at(upper_node).GetVars()(0);
            double x0dot = 0;
            double x1dot = 0;
            if (spline.at(lower_node).GetType() == FullDeriv) {
                x0dot = spline.at(lower_node).GetVars()(1);
            }
            if (spline.at(upper_node).GetType() == FullDeriv) {
                x1dot = spline.at(upper_node).GetVars()(1);
            }

            const double da2dt1 = -2*pow(deltat, -3)*(3*(x0 - x1) + deltat*(2*x0dot - x1dot)) -
                                  pow(deltat, -2)*(3*(x0 - x1) - (2*x0dot + x1dot));
            const double da3dt1 = 3*pow(deltat, -4)*(2*(x0 - x1) + deltat*(x0dot + x1dot)) +
                                  pow(deltat, -3)*(2*(x0 - x1) - (x0dot + x1dot));

            return da2dt1*pow(time_spline,2) + da3dt1* pow(time_spline, 3)/static_cast<double>(num_force_polys_);
        }

        return 0;
    }

    vector_t EndEffectorSplines::ComputeCoefPartialWrtTime(SplineType type, int coord, double time, int time_idx) const {
        const node_v& spline = SelectSpline(type, coord);
        const int upper_node = GetUpperNodeIdx(type, coord, time);
        const int lower_node = GetLowerNodeIdx(type, coord, time);

        double deltat = times_.at(upper_node).GetTime() - times_.at(lower_node).GetTime();
        if (deltat == 0) {
            deltat = times_.at(upper_node).GetTime() - times_.at(GetLowerNodeIdx(type, coord, time - 1e-4)).GetTime();
        }
        assert(deltat != 0);
        const double time_spline = time - times_.at(lower_node).GetTime();

        int vars_index, vars_affecting;
        std::tie(vars_index, vars_affecting) = GetVarsIdx(type, coord, time);

        vector_t coef_partials(vars_affecting);
        coef_partials.setZero();

        // check for direct dependence
        bool direct_dep = false;
        bool wrt_lower = false;

        if (type == Position) {
            direct_dep = true;

            const int node = ConvertContactNodeToSplineNode(time_idx);
            if (node == lower_node) {
                wrt_lower = true;
            }
        } else {
            const int node = ConvertContactNodeToSplineNode(time_idx);
            if (forces_.at(coord).at(node).GetType() == NoDeriv) {
                direct_dep = true;
            }

            if (direct_dep && node == lower_node) {
                wrt_lower = true;
            }
        }

        if (type == Force) {
            if (direct_dep) {
                if (spline.at(lower_node).GetType() == FullDeriv) {
                    assert(vars_affecting >= 2);
                    coef_partials(0) = Getx0CoefPartial(time_spline, deltat, wrt_lower);
                    coef_partials(1) = Getx0dotCoefPartial(time_spline, deltat, wrt_lower);
                    if (spline.at(upper_node).GetType() == FullDeriv) {
                        assert(vars_affecting == 4);
                        coef_partials(2) = Getx1CoefPartial(time_spline, deltat, !wrt_lower);
                        coef_partials(3) = Getx1dotCoefPartial(time_spline, deltat, !wrt_lower);
                    }
                } else {
                    assert(vars_affecting == 2);
                    coef_partials(0) = Getx1CoefPartial(time_spline, deltat, !wrt_lower);
                    coef_partials(1) = Getx1dotCoefPartial(time_spline, deltat, !wrt_lower);
                }

//                coef_partials(0) = Getx0CoefPartial(time, deltat, wrt_lower);
//                const int node = ConvertContactNodeToSplineNode(time_idx);
//                // Now need to determine if it is x0dot, or x1 next
//                if (spline.at(lower_node).GetType() == FullDeriv) {
//                    coef_partials(1) = Getx0dotCoefPartial(time, deltat, wrt_lower);
//                    coef_partials(2) = Getx1CoefPartial(time, deltat, !wrt_lower);
//                    if (spline.at(upper_node).GetType() == FullDeriv) {
//                        coef_partials(3) = Getx1dotCoefPartial(time, deltat, !wrt_lower);
//                    }
//                } else {
//                    coef_partials(1) = Getx1CoefPartial(time, deltat, !wrt_lower);
//                    if (spline.at(upper_node).GetType() == FullDeriv) {
//                        coef_partials(2) = Getx1dotCoefPartial(time, deltat, !wrt_lower);
//                    }
//                }
            } else {
                // TODO: Implement
            }
        } else {
            if (direct_dep) {
                if (spline.at(lower_node).GetType() == FullDeriv) {
                    assert(vars_affecting >= 2);
                    coef_partials(0) = Getx0CoefPartial(time_spline, deltat, wrt_lower);
                    coef_partials(1) = Getx0dotCoefPartial(time_spline, deltat, wrt_lower);
                    if (spline.at(upper_node).GetType() == NoDeriv) {
                        assert(vars_affecting == 4);
                        coef_partials(2) = Getx1CoefPartial(time_spline, deltat, !wrt_lower);
                        coef_partials(3) = Getx1dotCoefPartial(time_spline, deltat, !wrt_lower);
                    }
                } else if (vars_affecting == 1) {
                    coef_partials(0) = Getx0CoefPartial(time_spline, deltat, wrt_lower);
                } else {
                    assert(vars_affecting == 2);
                    coef_partials(0) = Getx0CoefPartial(time_spline, deltat, wrt_lower);
                    coef_partials(1) = Getx1CoefPartial(time_spline, deltat, wrt_lower);
                }
            } else {
                // TODO: Implement
            }
        }

        return coef_partials;
    }

    void EndEffectorSplines::SetVars(SplineType type, int coord, int node_idx, const vector_2t& vars) {
        // Force NoDerivs can only be zeros
        // Position NoDerivs need to be paired

        node_v& spline = SelectSpline(type, coord);
        if (spline.at(node_idx).GetType() == Empty) {
            throw std::runtime_error("Can't set this node's variables. This node is empty.");
        }

        if (type == Force && spline.at(node_idx).GetType() == NoDeriv) {
            throw std::runtime_error("Force spline cannot be changed at that node. Always set to 0.");
        }

        if (type == Position && coord != 2) {
            if ((node_idx < spline.size()-1 && forces_.at(coord).at(node_idx+1).GetType() == FullDeriv)) {
                spline.at(node_idx).SetVars(vars);
                spline.at(node_idx+spline_stride_).SetVars(vars);
            } else if(node_idx > 0 && forces_.at(coord).at(node_idx-1).GetType() == FullDeriv) {
                spline.at(node_idx).SetVars(vars);
                if (node_idx >= spline_stride_) {
                    spline.at(node_idx - spline_stride_).SetVars(vars);
                }
            } else {
                spline.at(node_idx).SetVars(vars);
            }
        } else if (type == Position) {
            if ((node_idx < spline.size()-1 && positions_.at(coord).at(node_idx).GetType() != FullDeriv &&
                    forces_.at(coord).at(node_idx+1).GetType() == FullDeriv)) {
                spline.at(node_idx).SetVars(vars);
                spline.at(node_idx+spline_stride_).SetVars(vars);
            } else if(node_idx > 0 && positions_.at(coord).at(node_idx).GetType() != FullDeriv &&
                    forces_.at(coord).at(node_idx-1).GetType() == FullDeriv) {
                spline.at(node_idx).SetVars(vars);
                if (node_idx >= spline_stride_) {
                    spline.at(node_idx - spline_stride_).SetVars(vars);
                }
            } else {
                spline.at(node_idx).SetVars(vars);
            }
        } else {
            spline.at(node_idx).SetVars(vars);
        }

    }

    void EndEffectorSplines::SetContactTimes(const time_v& contact_times) {

        assert(contact_times.size() == GetNumContacts()); // TODO: Check this

        const int num_nodes_start = GetNumNodes();
        int contact_idx = 0;
        for (int i = 0; i < GetNumNodes(); i++) {
            if (times_.at(i).GetType() == LiftOff || times_.at(i).GetType() == TouchDown) {
                assert(positions_.at(0).at(i).GetType() == NoDeriv);
                times_.at(i).SetTime(contact_times.at(contact_idx).GetTime());
                contact_idx++;
            } else if (forces_.at(0).at(i).GetType() == Empty) {
                times_.at(i).SetTime(times_.at(i-1).GetTime() + (contact_times.at(contact_idx).GetTime()
                                - contact_times.at(contact_idx-1).GetTime())/2);
            } else {
                double contact_time = 0.2 + contact_times.at(contact_idx - 1).GetTime();
                if (contact_idx < contact_times.size()) {
                    contact_time = contact_times.at(contact_idx).GetTime() - contact_times.at(contact_idx - 1).GetTime();
                }
                times_.at(i).SetTime(times_.at(i-1).GetTime() + contact_time/num_force_polys_);
            }
        }

        assert(num_nodes_start == GetNumNodes());
    }

    NodeType EndEffectorSplines::GetNodeType(SplineType type, int coord, int node_idx) const {
        const node_v& spline = SelectSpline(type, coord);
        return spline.at(node_idx).GetType();
    }

    int EndEffectorSplines::GetNumNodes() const {
        assert(forces_.at(0).size() == positions_.at(0).size());
        assert(forces_.at(0).size() == times_.size());
        return times_.size();
    }

    std::vector<int> EndEffectorSplines::GetMutableNodes(SplineType type, int coord) const {
        std::vector<int> mutable_nodes;
        switch (type) {
            case Force:
                for (int i = 0; i < GetNumNodes(); i++) {
                    if (forces_.at(coord).at(i).GetType() == FullDeriv) {
                        mutable_nodes.push_back(i);
                    }
                }
                return mutable_nodes;
            case Position:
                for (int i = 0; i < GetNumNodes(); i++) {
                    if (positions_.at(coord).at(i).GetType() != Empty) {
                        if (coord != 2) {
                            mutable_nodes.push_back(i);
                            if (i + spline_stride_ < GetNumNodes() &&
                                positions_.at(coord).at(i + spline_stride_).GetType() == NoDeriv) {
                                i += spline_stride_;
                            }
                        } else {
                            mutable_nodes.push_back(i);
//                            if (!(i - spline_stride_ >= 0 && positions_.at(coord).at(i - spline_stride_).GetType() == NoDeriv)) {
//
//                            } else
                            if (i + spline_stride_ < GetNumNodes() &&
                                    positions_.at(coord).at(i + spline_stride_).GetType() == NoDeriv) {
                                    i += spline_stride_;
                            }
                        }
                    }
                }
                return mutable_nodes;
            default:
                throw std::runtime_error("Unsupported spline type.");
        }
    }

    std::vector<double> EndEffectorSplines::GetTimes() const {
        std::vector<double> times;
        for (auto time : times_) {
            times.push_back(time.GetTime());
        }
        return times;
    }

    vector_t EndEffectorSplines::GetSplineAsQPVec(SplineType type, int coord) const {
        const node_v& spline = SelectSpline(type, coord);
        const std::vector<int> mut_nodes = GetMutableNodes(type, coord);

        int vec_size = 0;
        for (auto& it : mut_nodes) {
            if (spline.at(it).GetType() == NoDeriv) {
                vec_size++;
            } else {
                vec_size += 2;
            }
        }

        vector_t qp_vec(vec_size);
        qp_vec.setZero();
        int vec_idx = 0;
        for (auto& it : mut_nodes) {
            if (spline.at(it).GetType() == NoDeriv) {
                qp_vec(vec_idx) = spline.at(it).GetVars()(0);
                vec_idx++;
            } else {
                qp_vec.segment<2>(vec_idx) = spline.at(it).GetVars();
                vec_idx += 2;
            }
        }

        assert(vec_idx == vec_size);

        return qp_vec;
    }


    double EndEffectorSplines::GetEndTime() const {
        return times_.at(times_.size()-1).GetTime();
    }

    double EndEffectorSplines::GetStartTime() const {
        return times_.at(0).GetTime();
    }

    int EndEffectorSplines::GetTotalPolyVars(SplineType type, int coord) const {
        const node_v& spline = SelectSpline(type, 0);
        if (type == Force) {
            return 2*GetMutableNodes(type, coord).size();
        } else {
            return GetMutableNodes(type, coord).size();
        }
    }

    int EndEffectorSplines::GetNumContacts() const {
        int count = 0;
        for (int i = 0; i < GetNumNodes(); i++) {
            if (times_.at(i).GetType() == LiftOff || times_.at(i).GetType() == TouchDown) {
                count++;
            }
        }

        return count;
    }

    std::vector<double> EndEffectorSplines::GetContactTimeValues() const {
        std::vector<double> contact_times;
        for (int i = 0; i < GetNumNodes(); i++) {
            if (times_.at(i).GetType() == LiftOff || times_.at(i).GetType() == TouchDown) {
                assert(positions_.at(0).at(i).GetType() == NoDeriv);
                contact_times.push_back(times_.at(i).GetTime());
            }
        }

        return contact_times;
    }

    time_v EndEffectorSplines::GetContactTimes() const {
        time_v times;
        for (int i = 0; i < GetNumNodes(); i++) {
            if (times_.at(i).GetType() == LiftOff || times_.at(i).GetType() == TouchDown) {
                times.push_back(times_.at(i));
            }
        }

        return times;
    }

    int EndEffectorSplines::GetLowerNodeIdx(SplineType type, int coord, double time) const {
        if (time < times_.at(0).GetTime()) {
            throw std::runtime_error("Time requested is too small.");
        }

        if (time > times_.at(times_.size()-1).GetTime()) {
            throw std::runtime_error("Time requested is too large.");
        }

        const node_v& spline = SelectSpline(type, coord);

        for (int i = times_.size() - 1; i >= 0 ; i--) {
            if (time >= times_.at(i).GetTime() && spline.at(i).GetType() != Empty) {
                return i;
            }
        }

        throw std::runtime_error("Invalid time.");
    }

    int EndEffectorSplines::GetUpperNodeIdx(SplineType type, int coord, double time) const {
        if (time < times_.at(0).GetTime()) {
            throw std::runtime_error("Time requested is too small.");
        }

        if (time > times_.at(times_.size()-1).GetTime()) {
            throw std::runtime_error("Time requested is too large.");
        }

        const node_v& spline = SelectSpline(type, coord);

        for (int i = 0; i < times_.size(); i++) {
            if (time < times_.at(i).GetTime() && spline.at(i).GetType() != Empty) {
                return i;
            }
        }

        if (time == times_.at(times_.size()-1).GetTime()) {
            return times_.size()-1;
        }

        throw std::runtime_error("Invalid time.");
    }

    int EndEffectorSplines::ConvertContactNodeToSplineNode(int contact_idx) const {
        int contacts = 0;
        for (int i = 0; i < times_.size(); i++) {
            if (contacts == contact_idx) {
                return i;
            }
            if (forces_.at(0).at(i).GetType() == NoDeriv) {
                contacts++;
            }
        }

        return -1;
    }


    inline const EndEffectorSplines::node_v& EndEffectorSplines::SelectSpline(SplineType type, int coord) const {
        switch (type) {
            case Force:
                return forces_.at(coord);
            case Position:
                return positions_.at(coord);
            default:
                throw std::runtime_error("The provided spline type is not supported.");
                break;
        }
    }

    inline EndEffectorSplines::node_v& EndEffectorSplines::SelectSpline(SplineType type, int coord) {
        switch (type) {
            case Force:
                return forces_.at(coord);
            case Position:
                return positions_.at(coord);
            default:
                throw std::runtime_error("The provided spline type is not supported.");
                break;
        }
    }

    inline void EndEffectorSplines::UnsupportedSpline() {
        throw std::runtime_error("The provided spline type is not supported.");
    }

    double EndEffectorSplines::Getx0Coef(double time, double deltat) const {
        return 1 - (1 / pow(deltat, 2)) * 3 * pow(time, 2) +
               (1 / pow(deltat, 3)) * 2 * pow(time, 3);
    }

    double EndEffectorSplines::Getx1Coef(double time, double deltat) const {
        return (1 / pow(deltat, 2)) * 3 * pow(time, 2) -
               (1 / pow(deltat, 3)) * 2 * pow(time, 3);
    }

    double EndEffectorSplines::Getx0dotCoef(double time, double deltat) const {
        return time - (1 / deltat) * 2 * pow(time, 2) +
               (1 / pow(deltat, 2)) * pow(time, 3);
    }

    double EndEffectorSplines::Getx1dotCoef(double time, double deltat) const {
        return -(1 / deltat) * pow(time, 2) +
               (1 / pow(deltat, 2)) * pow(time, 3);
    }

    double EndEffectorSplines::Getx0CoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return -6*pow(DeltaT, -3)* pow(time, 2) + 6*pow(DeltaT, -2)*time
                   - 6*pow(DeltaT, 2)*pow(time,3) - 6*pow(DeltaT, 3)* pow(time, 2);
        } else {
            return 6*pow(DeltaT, -3)* pow(time, 2) + 6* pow(DeltaT, 2)* pow(time, 3);
        }
    }

    double EndEffectorSplines::Getx1CoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return 6*pow(DeltaT, -3)* pow(time, 2) - 6* pow(DeltaT, -2)*time
                   - 6* pow(DeltaT, -4)* pow(time, 3) + 6*pow(DeltaT, -3)* pow(time, 2);
        } else {
            return -6* pow(DeltaT, -3)* pow(time, 2) + 6* pow(DeltaT, -4)* pow(time, 3);
        }
    }

    double EndEffectorSplines::Getx0dotCoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return -1 -2*pow(DeltaT,-2)*pow(time, 2) + 4* pow(DeltaT, -1)*time
                   + 2* pow(DeltaT, -3)* pow(time, 3) - 3*pow(DeltaT, -2)*pow(time, 2);
        } else {
            return 2*pow(DeltaT, -2)* pow(time, 2) - 2* pow(DeltaT, -3)* pow(time, 3);
        }
    }

    double EndEffectorSplines::Getx1dotCoefPartial(double time, double DeltaT, bool wrt_t1) const {
        if (wrt_t1) {
            return -pow(DeltaT, -2)* pow(time, 2) + 2*pow(DeltaT, -1)*time
                   + 2*pow(DeltaT,-3)*pow(time, 3) - 3*pow(DeltaT, -2)*pow(time, 2);
        } else {
            return -pow(DeltaT, -2)* pow(time, 2) - 2*pow(DeltaT, -3)*pow(time, 3);
        }
    }

} // mpc