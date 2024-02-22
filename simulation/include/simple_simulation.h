//
// Created by zolkin on 2/21/24.
//

#ifndef BILEVEL_GAIT_GEN_SIMPLE_SIMULATION_H
#define BILEVEL_GAIT_GEN_SIMPLE_SIMULATION_H

#include <string>
#include <GLFW/glfw3.h>
#include <Eigen/Core>

#include "mujoco.h"

namespace simulation {
    // Simple single threaded simulation and visualizer
    class SimpleSimulation {
    public:
        SimpleSimulation(const std::string& robot_xml);

        mjtNum* GetControlPointer();

        mjData* GetDataPointer();

        const mjModel* GetModelPointer() const;

        void UpdateSim(double dt);

        void GetTrajViz(const std::vector<std::vector<Eigen::Vector3d>>& traj_viz,
                        const Eigen::Vector2d& box_sides,
                        const std::vector<Eigen::Vector2d>& box_centers);

    private:
        void UpdateTrajViz();

        std::vector<std::vector<Eigen::Vector3d>> traj_viz_;
        Eigen::Vector2d box_sides_;
        std::vector<Eigen::Vector2d> box_centers_;
    };
}


#endif //BILEVEL_GAIT_GEN_SIMPLE_SIMULATION_H
