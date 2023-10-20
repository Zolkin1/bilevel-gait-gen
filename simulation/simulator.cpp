//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <iostream>
#include "mujoco.h"

#include "simulator.h"

Simulator::Simulator() {
    std::cout << "Creating the simulator" << std::endl;

    // MuJoCo data structures
    mjModel* m = NULL;                  // MuJoCo model
    mjData* d = NULL;                   // MuJoCo data
    mjvCamera cam;                      // abstract camera
    mjvOption opt;                      // visualization options
    mjvScene scn;                       // abstract scene
    mjrContext con;                     // custom GPU context
}

void Simulator::SetupSimulator() {

}