//
// Created by zolkin on 1/30/24.
//

#ifndef BILEVEL_GAIT_GEN_END_EFFECTOR_SPLINES_H
#define BILEVEL_GAIT_GEN_END_EFFECTOR_SPLINES_H

#include "spline_node.h"

namespace mpc {
    class EndEffectorSplines {
        using node_v = std::vector<SplineNode>;
        using time_v = std::vector<double>;
        using vector_4t = Eigen::Vector4d;
    public:
        enum SplineType {
            Force = 0,
            Position = 1
        };

        EndEffectorSplines(int num_contacts,
                           const std::vector<double>& times,
                           bool start_on_constant,
                           int num_force_polys);

        double ValueAt(SplineType type, double time) const;

        // ------------------ Setters ------------------ //
        void SetVars(SplineType type, int node_idx, const vector_2t& vars);


        // ------------------ Getters ------------------ //
        NodeType GetNodeType(SplineType type, int node_idx) const;

        int GetNumNodes() const;        // Nodes (includes intermediate nodes)

        // Function list:
        // - Spline Value
        // - Vars coefficients -- Compressed view
        // - Vars index -- Compressed view
        // - Num vars effecting -- Compressed view
        // - Setting vars
        // - Total vars
        // - Num nodes in the spline
        // - Num of each type

    protected:
    private:
        int GetLowerNodeIdx(SplineType type, double time) const;
        int GetUpperNodeIdx(SplineType type, double time) const;

        node_v& SelectSpline(SplineType type);

        const node_v& SelectSpline(SplineType type) const;


        static void UnsupportedSpline() ;

        node_v forces_;
        node_v positions_;
        time_v times_;
        const int num_force_polys_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_END_EFFECTOR_SPLINES_H
