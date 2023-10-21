//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SIMULATOR_H
#define BILEVEL_GAIT_GEN_SIMULATOR_H

#include <string>
#include <memory>

#include "mujoco.h"
#include "GLFW/glfw3.h"

#include "simulation_robot.h"

namespace simulator {
    // MuJoCo data structures
    extern mjModel* model_;                  // MuJoCo model
    extern mjData* data_;                    // MuJoCo data
    extern mjvCamera cam_;                   // abstract camera
    extern mjvOption opt_;                   // visualization options
    extern mjvScene scene_;                  // abstract scene
    extern mjrContext context_;              // custom GPU context


    // mouse interaction
    extern bool button_left_;
    extern bool button_middle_;
    extern bool button_right_;
    extern double lastx_;
    extern double lasty_;

    /**
     * Simulator class.
     * Should be able to take in a robot model and a controller and simulate it.
     * Simulator will handle all the Mujoco interfaces itself.
     *
     * Simulator will need access to:
     * - Robot visualization information (like meshes potentially)
     * - Controllers
     * - Robot dynamics (like a urdf that can be passed to Mujoco)
     *
     * This simulator is purely for verification purposes as the MPC will have an internal simulator.
     * Therefore the simulator will always want to open up a seperate visualziation window.
     */
    class Simulator {
    public:
        /**
         * Simulator constructor.
         */
        Simulator(std::unique_ptr<SimulationRobot>& robot, double viz_rate);

        /**
         * Setup the simulator
         */
        void SetupSimulator();

        void RunSimulator();

        /**
         * Interfaces with the mujoco data to set the state of the robot.
         * Called in SetupSimulator.
         */
        void SetState(const Eigen::VectorXd& config, const Eigen::VectorXd& velocities);

    private:
        // mouse move callback
        static void MouseMove(GLFWwindow *window, double xpos, double ypos);
        static void MouseButton(GLFWwindow* window, int button, int act, int mods);
        static void Scroll(GLFWwindow* window, double xoffset, double yoffset);

        GLFWwindow* window_;

        double visualization_rate_;

        std::unique_ptr<SimulationRobot> robot_;
    };
}   // simulator


#endif //BILEVEL_GAIT_GEN_SIMULATOR_H
