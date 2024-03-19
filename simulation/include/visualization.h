//
// Created by zolkin on 12/12/23.
//

#ifndef BILEVEL_GAIT_GEN_VISUALIZATION_H
#define BILEVEL_GAIT_GEN_VISUALIZATION_H

#include <Eigen/Core>
#include <GLFW/glfw3.h>

#include "mujoco.h"

namespace simulation {
    using vector_t = Eigen::VectorXd;
    using matrix_t =  Eigen::MatrixXd;

    class Visualizer {
    public:
        Visualizer(const std::string& robot_xml);

        void UpdateState(const vector_t& config);
        void UpdateViz(double dt);
        void GetTrajViz(const std::vector<std::vector<Eigen::Vector3d>>& traj_viz,
                        const Eigen::Vector2d& box_sides,
                        const std::vector<Eigen::Vector2d>& box_centers);

        const mjModel* GetModel() const;
    protected:
    private:
        void UpdateTrajViz();

        std::vector<std::vector<Eigen::Vector3d>> traj_viz_;
        Eigen::Vector2d box_sides_;
        std::vector<Eigen::Vector2d> box_centers_;
    };
} // simulation


#endif //BILEVEL_GAIT_GEN_VISUALIZATION_H
