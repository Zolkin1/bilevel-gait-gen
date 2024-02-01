//
// Created by zolkin on 1/30/24.
//

#include "spline_node.h"

namespace mpc {
    SplineNode::SplineNode(mpc::NodeType type, int node_idx, const vector_2t& vars)
    : type_(type), node_idx_(node_idx), vars_(vars) {
        if (type_ == NoDeriv && vars_(1) != 0) {
            throw std::runtime_error("Type does not match the provided variables. Expecting no derivative term.");
        }
    }

    NodeType SplineNode::GetType() const {
        return type_;
    }

    vector_2t SplineNode::GetVars() const {
        if (type_ == Empty) {
            throw std::runtime_error("Can't get the vars in this node. This node is set to empty.");
        }

        return vars_;
    }

    int SplineNode::GetNodeIdx() const {
        return node_idx_;
    }

    void SplineNode::SetVars(const vector_2t& vars) {
        if (type_ == Empty) {
            throw std::runtime_error("Can't set the vars in this node. This node is set to empty.");
        }

        if (type_ == NoDeriv) {
            vars_(0) = vars(0);
        } else {
            vars_ = vars;
        }
    }

    SplineNode::SplineNode(const SplineNode& node)
    : type_(node.type_), node_idx_(node.node_idx_), vars_(node.vars_) {}

    SplineNode& SplineNode::operator=(const SplineNode& node) {
        if (this == &node) {
            return *this;
        }

        vars_ = node.vars_;
        type_ = node.type_;
        node_idx_ = node.node_idx_;
        return *this;
    }

} // mpc