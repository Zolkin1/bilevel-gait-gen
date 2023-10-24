//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SIMULATOR_H
#define BILEVEL_GAIT_GEN_SIMULATOR_H

#include <string>
#include <memory>
#include <mutex>
#include <optional>
#include <cstdarg>

#include "simulate.h"
#include "mujoco.h"
#include "GLFW/glfw3.h"

#include "simulation_robot.h"

namespace simulator {
    // MuJoCo data structures
    extern mjModel* model_;                  // MuJoCo model
    extern mjData* data_;                    // MuJoCo data
    extern mjvCamera cam_;                   // abstract camera
    extern mjvOption opt_;                   // visualization options
    extern mjvPerturb pert;


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
     *
     * Simulator will make one thread for rendering and one thread for simulating.
     */
    class Simulator {
    public:
        /**
         * Simulator constructor.
         */
        explicit Simulator(std::unique_ptr<SimulationRobot>& robot);

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

        static void SimulateLoop(mujoco::Simulate* sim, SimulationRobot* robot);

        // ----------------------- Member Variables ----------------------- //
        std::unique_ptr<SimulationRobot> robot_;

        std::unique_ptr<mujoco::Simulate> sim;

        static constexpr double syncMisalign = 0.1;        // maximum mis-alignment before re-sync (simulation seconds)
        static constexpr double simRefreshFraction = 0.7;  // fraction of refresh available for simulation


    }; // Simulator

}   // simulator


#endif //BILEVEL_GAIT_GEN_SIMULATOR_H
