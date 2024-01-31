//
// Created by zolkin on 1/30/24.
//
#include <assert.h>

#include "end_effector_splines.h"

namespace mpc {

    EndEffectorSplines::EndEffectorSplines(int num_contacts,
                                           const std::vector<double>& times,
                                           bool start_on_constant,
                                           int num_force_polys)
                                           : num_force_polys_(num_force_polys) {
        const vector_2t ZERO_VARS = vector_2t::Zero();

        std::vector<NodeType> force_type_pattern;
        std::vector<NodeType> position_type_pattern;
        for (int i = 0; i < num_force_polys_ + 1; i++) {
            if (start_on_constant) {
                if (i < 2) {
                    force_type_pattern.push_back(NoDeriv);
                    position_type_pattern.push_back(NoDeriv);
                } else {
                    force_type_pattern.push_back(FullDeriv);
                    position_type_pattern.push_back(Empty);
                }
            } else {
                if (i < num_force_polys) {
                    force_type_pattern.push_back(NoDeriv);
                    position_type_pattern.push_back(NoDeriv);
                } else {
                    force_type_pattern.push_back(FullDeriv);
                    position_type_pattern.push_back(Empty);
                }
            }
        }

        int i = 0;
        int j = 0;
        while (i < num_contacts) {
            forces_.emplace_back(force_type_pattern.at(j%force_type_pattern.size()), j, ZERO_VARS);
            positions_.emplace_back(position_type_pattern.at(j%position_type_pattern.size()), j, ZERO_VARS);

            if (force_type_pattern.at(j%position_type_pattern.size()) == FullDeriv) {
                times_.push_back(times_.at(times_.size()-1) + (times.at(i) - times.at(i - 1))/num_force_polys);
            } else {
                times_.push_back(times.at(i));
                i++;
            }
            j++;
        }

        assert(forces_.size() == positions_.size());
        assert(forces_.size() == times_.size());
    }

    double EndEffectorSplines::ValueAt(SplineType type, double time) const {
        const node_v& spline = SelectSpline(type);
        const int lower_node = GetLowerNodeIdx(type, time);
        const int upper_node = GetUpperNodeIdx(type, time);

        const double deltat = times_.at(upper_node) - times_.at(lower_node);
        const double time_spline = time - times_.at(lower_node);
        const vector_4t vars {spline.at(lower_node).GetVars()(0),         // x0
                              spline.at(upper_node).GetVars()(0),       // x1
                              spline.at(lower_node).GetVars()(1),         // x0dot
                              spline.at(upper_node).GetVars()(1)};      // x1dot

        const double a2 = -(1 / pow(deltat, 2)) * 3 * (vars(0) - vars(1)) -
                    (1 / deltat) * (2 * vars(2) + vars(3));
        const double a3 = (1 / pow(deltat, 3)) * 2 * (vars(0) - vars(1)) +
                    (1 / pow(deltat, 2)) * (vars(2) + vars(3));

        const double val = vars(0) + vars(2) * time_spline + a2 * pow(time_spline, 2) + a3 * pow(time_spline, 3);

        return val;

    }

    void EndEffectorSplines::SetVars(SplineType type, int node_idx, const vector_2t& vars) {
        node_v& spline = SelectSpline(type);
        spline.at(node_idx).SetVars(vars);
    }

    NodeType EndEffectorSplines::GetNodeType(SplineType type, int node_idx) const {
        const node_v& spline = SelectSpline(type);
        return spline.at(node_idx).GetType();
    }

    int EndEffectorSplines::GetNumNodes() const {
        assert(forces_.size() == positions_.size());
        assert(forces_.size() == times_.size());
        return times_.size();
    }

    int EndEffectorSplines::GetLowerNodeIdx(SplineType type, double time) const {
        if (time < times_.at(0)) {
            throw std::runtime_error("Time requested is too small.");
        }

        if (time > times_.at(times_.size()-1)) {
            throw std::runtime_error("Time requested is too large.");
        }

        const node_v& spline = SelectSpline(type);

        for (int i = times_.size() - 1; i >= 0 ; i--) {
            if (time >= times_.at(i) && spline.at(i).GetType() != Empty) {
                return i;
            }
        }

        throw std::runtime_error("Invalid time.");
    }

    int EndEffectorSplines::GetUpperNodeIdx(SplineType type, double time) const {
        if (time < times_.at(0)) {
            throw std::runtime_error("Time requested is too small.");
        }

        if (time >= times_.at(times_.size()-1)) {
            throw std::runtime_error("Time requested is too large.");
        }

        const node_v& spline = SelectSpline(type);

        for (int i = 0; i < times_.size(); i++) {
            if (time < times_.at(i) && spline.at(i).GetType() != Empty) {
                return i;
            }
        }

        throw std::runtime_error("Invalid time.");
    }


    inline const EndEffectorSplines::node_v& EndEffectorSplines::SelectSpline(mpc::EndEffectorSplines::SplineType type) const {
        switch (type) {
            case Force:
                return forces_;
            case Position:
                return positions_;
            default:
                UnsupportedSpline();
                break;
        }
    }

    inline EndEffectorSplines::node_v& EndEffectorSplines::SelectSpline(mpc::EndEffectorSplines::SplineType type) {
        switch (type) {
            case Force:
                return forces_;
            case Position:
                return positions_;
            default:
                UnsupportedSpline();
                break;
        }
    }

    inline void EndEffectorSplines::UnsupportedSpline() {
        throw std::runtime_error("The provided spline type is not supported.");
    }

} // mpc