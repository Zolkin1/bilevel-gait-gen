//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

/*
 * Demo file. This should show how to setup the controller and connect it to a simulation.
 * I should be able to take in a urdf file (or maybe a mujoco file) and create the model, create the controller,
 * create the simulation and run it all.
 */

#include <iostream>
#include "simulator.h"

int main() {
    std::cout << "Hello world!" << std::endl;

    Simulator sim;

    sim.SetupSimulator();

    // Read in robot file (command line argument)

    // Create a model to be used by MPC (this will let the user specify anything not in the file)
    // this class will wrap pinocchio and provide everything the MPC needs

    // Create the controller
    // - Define cost functions and constraints
    // - Specify how the constraints will be treated

    // Initialize the simulation

    // Run the controller and retrieve diagnostics

    return 0;
}
