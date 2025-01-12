//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <iostream>
#include <thread>
#include <chrono>
#include <mutex>

#include "mujoco.h"
#include "simulate.h"
#include "glfw_adapter.h"

#include "simulator.h"
namespace simulator {
    // MuJoCo data structures
    mjModel* model_ = nullptr;                 // MuJoCo model
    mjData* data_ = nullptr;                   // MuJoCo data
    mjvCamera cam_;                         // abstract camera
    mjvOption opt_;                         // visualization options
    mjvPerturb pert;

    using Seconds = std::chrono::duration<double>;

    Simulator::Simulator(std::unique_ptr<SimulationRobot>& robot) {
        char error[1000] = "Could not load binary model";
        model_ = mj_loadXML(robot->GetRobotXMLFile().c_str(), nullptr, error, 1000);

        if (!model_) {
            mju_error("Load model error: %s", error);
        }

        robot->SetSimModel(model_);
    }

    void Simulator::SetupSimulator(const std::unique_ptr<SimulationRobot>& robot) {
        // initialize visualization data structures
        mjv_defaultCamera(&cam_);
        mjv_defaultOption(&opt_);
        mjv_defaultPerturb(&pert);

        // create a mujoco simulate object
        sim = std::make_unique<mujoco::Simulate>(
                std::make_unique<mujoco::GlfwAdapter>(),
                &cam_, &opt_, &pert, /* is_passive = */ false
        );

        sim->LoadMessage(robot->GetRobotXMLFile().c_str());
        sim->user_scn = new mjvScene;
        mjv_defaultScene(sim->user_scn);
        mjv_makeScene(model_, sim->user_scn, 2000);

        data_ = mj_makeData(model_);

        // Set state to current initial condition
        this->SetState(robot->ConvertPinocchioConfigToMujoco(robot->GetInitConfig()),
                       robot->ConvertPinocchioVelToMujoco(robot->GetInitVelocities()));

    }

    void Simulator::RunSimulator(std::unique_ptr<SimulationRobot>& robot) {

        // Create a simulation thread
        std::thread simulationThreadHandle(&SimulateLoop, sim.get(), robot.get());

        sim->RenderLoop();
        simulationThreadHandle.join();
    }

    void Simulator::SetState(const Eigen::VectorXd& config, const Eigen::VectorXd& velocities) {
        if (config.size() != model_->nq) {
            std::cerr << "Trying to initialize the Mujoco configuration with the wrong sized vector." <<
                        " Need to a vector of size " << model_->nq << ", the vector provided is of size " <<
                        config.size() <<"." << "Ignoring this set state." << std::endl;
            return;
        }

        if (velocities.size() != model_->nv) {
            std::cerr << "Trying to initialize the Mujoco velocities with the wrong sized vector." <<
                      " Need to a vector of size " << model_->nv << ", the vector provided is of size " <<
                      velocities.size() <<"." << "Ignoring this set state." << std::endl;
            return;
        }

        for (int i = 0; i < model_->nq; i++) {
            data_->qpos[i] = config(i);
        }

        for (int i = 0; i < model_->nv; i++) {
            data_->qvel[i] = velocities(i);
        }

    }

    void Simulator::SimulateLoop(mujoco::Simulate* sim, SimulationRobot* robot) {
        sim->Load(model_, data_, robot->GetRobotXMLFile().c_str());
        mj_forward(model_, data_);

        // cpu-sim syncronization point
        std::chrono::time_point<mujoco::Simulate::Clock> syncCPU;
        mjtNum syncSim = 0;

        mjtNum last_control_time = data_->time;

        // run until asked to exit
        while (!sim->exitrequest.load()) {
            // sleep for 1 ms or yield, to let main thread run
            //  yield results in busy wait - which has better timing but kills battery life
            if (sim->run && sim->busywait) {
                std::this_thread::yield();
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }

            {
                // lock the sim mutex
                const std::unique_lock<std::recursive_mutex> lock(sim->mtx);

                // run only if model is present
                if (model_) {
                    // running
                    if (sim->run) {
                        bool stepped = false;

                        // record cpu time at start of iteration
                        const auto startCPU = mujoco::Simulate::Clock::now();

                        // elapsed CPU and simulation time since last sync
                        const auto elapsedCPU = startCPU - syncCPU;
                        double elapsedSim = data_->time - syncSim;

                        // requested slow-down factor
                        double slowdown = 100 / mujoco::Simulate::percentRealTime[sim->real_time_index];

                        // misalignment condition: distance from target sim time is bigger than syncmisalign
                        bool misaligned =
                                mju_abs(Seconds(elapsedCPU).count()/slowdown - elapsedSim) > syncMisalign;

                        // out-of-sync (for any reason): reset sync times, step
                        if (elapsedSim < 0 || elapsedCPU.count() < 0 || syncCPU.time_since_epoch().count() == 0 ||
                            misaligned || sim->speed_changed) {
                            // re-sync
                            syncCPU = startCPU;
                            syncSim = data_->time;
                            sim->speed_changed = false;

                            // Get control input - if the time is less than the last time, a reset happened, so apply a control
                            if (data_->time - last_control_time >= robot->GetController()->GetRate() ||
                                data_->time < last_control_time) {

                                robot->GetControlAction(data_, data_->ctrl);
                                std::vector<std::vector<Eigen::Vector3d>> viz_data = robot->GetTrajViz();
                                UpdateVizGeoms(sim, viz_data);

                                last_control_time = data_->time;
                            }
                            // run single step, let next iteration deal with timing
                            mj_step(model_, data_);
                            stepped = true;
                        } else { // in-sync: step until ahead of cpu
                            bool measured = false;
                            mjtNum prevSim = data_->time;

                            double refreshTime = simRefreshFraction / sim->refresh_rate;

                            // step while sim lags behind cpu and within refreshTime
                            while (Seconds((data_->time - syncSim)*slowdown) < mujoco::Simulate::Clock::now() - syncCPU &&
                                   mujoco::Simulate::Clock::now() - startCPU < Seconds(refreshTime)) {
                                // measure slowdown before first step
                                if (!measured && elapsedSim) {
                                    sim->measured_slowdown =
                                            std::chrono::duration<double>(elapsedCPU).count() / elapsedSim;
                                    measured = true;
                                }

                                // Get control input - if the time is less than the last time, a reset happened, so apply a control
                                if (data_->time - last_control_time >= robot->GetController()->GetRate() ||
                                    data_->time < last_control_time) {

                                    robot->GetControlAction(data_, data_->ctrl);
                                    std::vector<std::vector<Eigen::Vector3d>> viz_data = robot->GetTrajViz();
                                    UpdateVizGeoms(sim, viz_data);

                                    last_control_time = data_->time;
                                }
                                // call mj_step
                                mj_step(model_, data_);
                                stepped = true;

                                // break if reset
                                if (data_->time < prevSim) {
                                    break;
                                }
                            }
                        }

                        // save current state to history buffer
                        if (stepped) {
                            sim->AddToHistory();
                        }
                    }

                        // paused
                    else {
                        // run mj_forward, to update rendering and joint sliders
                        mj_forward(model_, data_);
                        sim->speed_changed = true;
                    }
                }
            }  // release std::lock_guard<std::mutex>
        }
    }

    void Simulator::UpdateVizGeoms(mujoco::Simulate* sim, const std::vector<std::vector<Eigen::Vector3d>>& viz_data) {
        if (!viz_data.empty()) {
            sim->user_scn->ngeom = 0;
            // TODO: Free other geoms?

            const mjtNum mat[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
            const float rgba_va[4] = {1, 1, 0, 1};
            const mjtNum size[3] = {1, 1, 1};
            const mjtNum pos[3] = {1, 1, 1};
            mjtNum from[3];
            mjtNum to[3];

            for (int ee = 0; ee < viz_data.size(); ee++) { // end effectors
                for (int node = 0; node < viz_data.at(ee).size() - 1; node++) {   // node
                    for (int coord = 0; coord < 3; coord++) {   // coordinate
                        from[coord] = viz_data.at(ee).at(node)(coord);
                        to[coord] = viz_data.at(ee).at(node+1)(coord);
                    }

                    sim->user_scn->ngeom += 1;
                    mjv_initGeom(&sim->user_scn->geoms[sim->user_scn->ngeom - 1],
                                 mjtGeom::mjGEOM_CAPSULE,
                                 size,
                                 pos,
                                 nullptr,
                                 rgba_va);

                    mjv_connector(&sim->user_scn->geoms[sim->user_scn->ngeom - 1],
                                  mjtGeom::mjGEOM_CAPSULE,
                                  0.005,
                                  from,
                                  to);
                }
            }
        }
    }

} // simulator