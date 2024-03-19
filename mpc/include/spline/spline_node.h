//
// Created by zolkin on 1/30/24.
//

#ifndef BILEVEL_GAIT_GEN_SPLINE_NODE_H
#define BILEVEL_GAIT_GEN_SPLINE_NODE_H

#include <Eigen/Core>

namespace mpc {
    using vector_t = Eigen::VectorXd;
    using vector_2t = Eigen::Vector2d;

    enum NodeType {
        NoDeriv = 0,
        FullDeriv = 1,
        Empty = 2
    };
    class SplineNode {
    public:
        SplineNode(NodeType type, int node_idx, const vector_2t& vars);

        SplineNode(const SplineNode& node);

        SplineNode& operator=(const SplineNode& node);

        NodeType GetType() const;
        vector_2t GetVars() const;
        int GetNodeIdx() const;

        void SetVars(const vector_2t& vars);

    protected:
    private:
        NodeType type_;   // Node type
        int node_idx_;    // Where I am in the spline
        vector_2t vars_;        // position and derivative that describe the spline at the point
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_SPLINE_NODE_H
