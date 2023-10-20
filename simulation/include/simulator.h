//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SIMULATOR_H
#define BILEVEL_GAIT_GEN_SIMULATOR_H

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
     // TODO: Take in robot model
    Simulator();

    /**
     * Setup the simulator
     */
    void SetupSimulator();
};


#endif //BILEVEL_GAIT_GEN_SIMULATOR_H
