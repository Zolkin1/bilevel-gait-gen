//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <iostream>
#include <functional>

#include "mujoco.h"
#include "GLFW/glfw3.h"

#include "simulator.h"
namespace simulator {
    // MuJoCo data structures
    mjModel* model_ = NULL;                 // MuJoCo model
    mjData* data_ = NULL;                   // MuJoCo data
    mjvCamera cam_;                         // abstract camera
    mjvOption opt_;                         // visualization options
    mjvScene scene_;                        // abstract scene
    mjrContext context_;                    // custom GPU context


    // mouse interaction
    bool button_left_ = false;
    bool button_middle_ = false;
    bool button_right_ = false;
    double lastx_ = 0.0;
    double lasty_ = 0.0;

    Simulator::Simulator(std::unique_ptr<SimulationRobot>& robot, double viz_rate) : robot_(std::move(robot)) {
        char error[1000] = "Could not load binary model";
        model_ = mj_loadXML(robot_->GetRobotXMLFile().c_str(), 0, error, 1000);

        if (!model_) {
            mju_error("Load model error: %s", error);
        }

        // make data
        data_ = mj_makeData(model_);

        if (!glfwInit()) {
            mju_error("Could not initialize GLFW");
        }

        // create window, make OpenGL context current, request v-sync
        window_ = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
        glfwMakeContextCurrent(window_);
        glfwSwapInterval(1);

        visualization_rate_ = viz_rate;
    }

    void Simulator::SetupSimulator() {
        // initialize visualization data structures
        mjv_defaultCamera(&cam_);
        mjv_defaultOption(&opt_);
        mjv_defaultScene(&scene_);
        mjr_defaultContext(&context_);

        // create scene and context
        mjv_makeScene(model_, &scene_, 2000);
        mjr_makeContext(model_, &context_, mjFONTSCALE_150);

        glfwSetCursorPosCallback(window_, Simulator::MouseMove);
        glfwSetMouseButtonCallback(window_, MouseButton);
        glfwSetScrollCallback(window_, Scroll);

        // Set state to current initial condition
        this->SetState(robot_->GetConfig(), robot_->GetVelocities());

        // Set the control callback function
        //auto ctrl = [capture0 = robot_.get()](auto && PH1, auto && PH2) { capture0->GetControlAction(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2)); };
        //mjcb_control = ctrl;
    }

    void Simulator::RunSimulator() {
        while (!glfwWindowShouldClose(window_)) {
            // advance interactive simulation for 1/60 sec
            //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
            //  this loop will finish on time for the next frame to be rendered at 60 fps.
            //  Otherwise add a cpu timer and exit this loop when it is time to render.
            mjtNum simstart = data_->time;
            while (data_->time - simstart < visualization_rate_) {
                mj_step(model_, data_);
            }

            // get framebuffer viewport
            mjrRect viewport = {0, 0, 0, 0};
            glfwGetFramebufferSize(window_, &viewport.width, &viewport.height);

            // update scene and render
            mjv_updateScene(model_, data_, &opt_, NULL, &cam_, mjCAT_ALL, &scene_);
            mjr_render(viewport, &scene_, &context_);

            // swap OpenGL buffers (blocking call due to v-sync)
            glfwSwapBuffers(window_);

            // process pending GUI events, call GLFW callbacks
            glfwPollEvents();
        }

        //free visualization storage
        mjv_freeScene(&scene_);
        mjr_freeContext(&context_);

        // free MuJoCo model and data
        mj_deleteData(data_);
        mj_deleteModel(model_);

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

    // ----------------------- GLFW Interaction Functions ----------------------- //
    void Simulator::MouseMove(GLFWwindow *window, double xpos, double ypos) {
        // no buttons down: nothing to do
        if (!button_left_ && !button_middle_ && !button_right_) {
            return;
        }

        // compute mouse displacement, save
        double dx = xpos - lastx_;
        double dy = ypos - lasty_;
        lastx_ = xpos;
        lasty_ = ypos;

        // get current window size
        int width, height;
        glfwGetWindowSize(window, &width, &height);

        // get shift key state
        bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                          glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

        // determine action based on mouse button
        mjtMouse action;
        if (button_right_) {
            action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
        } else if (button_left_) {
            action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
        } else {
            action = mjMOUSE_ZOOM;
        }

        // move camera
        mjv_moveCamera(model_, action, dx / height, dy / height, &scene_, &cam_);
    }

    // mouse button callback
    void Simulator::MouseButton(GLFWwindow* window, int button, int act, int mods) {
        // update button state
        button_left_ = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
        button_middle_ = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
        button_right_ = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

        // update mouse position
        glfwGetCursorPos(window, &lastx_, &lasty_);
    }

    // scroll callback
    void Simulator::Scroll(GLFWwindow* window, double xoffset, double yoffset) {
        // emulate vertical mouse motion = 5% of window height
        mjv_moveCamera(model_, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scene_, &cam_);
    }
}