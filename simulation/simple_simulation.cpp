//
// Created by zolkin on 2/21/24.
//

#include "simple_simulation.h"
#include "timer.h"

#include <cstdio>
#include <cstring>

#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>
#include <chrono>

namespace simulation {

// MuJoCo data structures
    mjModel *m = NULL;                  // MuJoCo model
    mjData *d = NULL;                   // MuJoCo data
    mjvCamera cam;                      // abstract camera
    mjvOption opt;                      // visualization options
    mjvScene scn;                       // abstract scene
    mjrContext con;                     // custom GPU context

    GLFWwindow* window;

// mouse interaction
    bool button_left = false;
    bool button_middle = false;
    bool button_right = false;
    double lastx = 0;
    double lasty = 0;

    bool pause = false;


// keyboard callback
    void keyboard(GLFWwindow *window, int key, int scancode, int act, int mods) {
        // backspace: reset simulation
        if (act == GLFW_PRESS && key == GLFW_KEY_BACKSPACE) {
            mj_resetData(m, d);
            mj_forward(m, d);
        }

        if (act == GLFW_PRESS && key == GLFW_KEY_SPACE) {
            pause = !pause;
        }
    }


// mouse button callback
    void mouse_button(GLFWwindow *window, int button, int act, int mods) {
        // update button state
        button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
        button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS);
        button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);

        // update mouse position
        glfwGetCursorPos(window, &lastx, &lasty);
    }


// mouse move callback
    void mouse_move(GLFWwindow *window, double xpos, double ypos) {
        // no buttons down: nothing to do
        if (!button_left && !button_middle && !button_right) {
            return;
        }

        // compute mouse displacement, save
        double dx = xpos - lastx;
        double dy = ypos - lasty;
        lastx = xpos;
        lasty = ypos;

        // get current window size
        int width, height;
        glfwGetWindowSize(window, &width, &height);

        // get shift key state
        bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                          glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

        // determine action based on mouse button
        mjtMouse action;
        if (button_right) {
            action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
        } else if (button_left) {
            action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
        } else {
            action = mjMOUSE_ZOOM;
        }

        // move camera
        mjv_moveCamera(m, action, dx / height, dy / height, &scn, &cam);
    }


// scroll callback
    void scroll(GLFWwindow *window, double xoffset, double yoffset) {
        // emulate vertical mouse motion = 5% of window height
        mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05 * yoffset, &scn, &cam);
    }


// main function
    int main(int argc, const char **argv) {
        // check command-line arguments
        if (argc != 2) {
            std::printf(" USAGE:  basic modelfile\n");
            return 0;
        }



        // run main loop, target real-time simulation and 60 fps rendering
        while (!glfwWindowShouldClose(window)) {
            // advance interactive simulation for 1/60 sec
            //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
            //  this loop will finish on time for the next frame to be rendered at 60 fps.
            //  Otherwise add a cpu timer and exit this loop when it is time to render.
            mjtNum simstart = d->time;
            while (d->time - simstart < 1.0 / 60.0) {
                mj_step(m, d);
            }

            // get framebuffer viewport
            mjrRect viewport = {0, 0, 0, 0};
            glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

            // update scene and render
            mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
            mjr_render(viewport, &scn, &con);

            // swap OpenGL buffers (blocking call due to v-sync)
            glfwSwapBuffers(window);

            // process pending GUI events, call GLFW callbacks
            glfwPollEvents();
        }

        //free visualization storage
        mjv_freeScene(&scn);
        mjr_freeContext(&con);

        // free MuJoCo model and data
        mj_deleteData(d);
        mj_deleteModel(m);

        // terminate GLFW (crashes with Linux NVidia drivers)
#if defined(__APPLE__) || defined(_WIN32)
        glfwTerminate();
#endif

        return 1;
    }

    SimpleSimulation::SimpleSimulation(const std::string& robot_xml) {
        // load and compile model
        char error[1000] = "Could not load binary model";
        m = mj_loadXML(robot_xml.c_str(), 0, error, 1000);

        if (!m) {
            mju_error("Load model error: %s", error);
        }

        // make data
        d = mj_makeData(m);

        // init GLFW
        if (!glfwInit()) {
            mju_error("Could not initialize GLFW");
        }

        // create window, make OpenGL context current, request v-sync
        window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);

        // initialize visualization data structures
        mjv_defaultCamera(&cam);
        mjv_defaultOption(&opt);
        mjv_defaultScene(&scn);
        mjr_defaultContext(&con);

        // create scene and context
        mjv_makeScene(m, &scn, 2000);
        mjr_makeContext(m, &con, mjFONTSCALE_150);

        // install GLFW mouse and keyboard callbacks
        glfwSetKeyCallback(window, keyboard);
        glfwSetCursorPosCallback(window, mouse_move);
        glfwSetMouseButtonCallback(window, mouse_button);
        glfwSetScrollCallback(window, scroll);

        run_num_ = 0;
    }

    void SimpleSimulation::UpdateSim(double dt) {
        using namespace std::chrono;
        steady_clock::time_point t1 = steady_clock::now();

        duration<double> time_span(0);
        double time = d->time;
        while(!glfwWindowShouldClose(window) && time_span.count() < dt) {

            if (!pause) {
                // Allows for multiple steps before the controller is re-computed

                // Set for debugging


                utils::Timer timer("sim physics");
                timer.StartTimer();
                while (d->time < time + dt) {
//                    d->qpos[0] = 0;
//                    d->qpos[1] = 0;
//                    d->qpos[2] = 0.34;
//                    d->qpos[3] = 1;
//                    d->qpos[4] = 0;
//                    d->qpos[5] = 0;
//                    d->qpos[6] = 0;
//
//                    d->qvel[0] = 0;
//                    d->qvel[1] = 0;
//                    d->qvel[2] = 0;
//                    d->qvel[3] = 0;
//                    d->qvel[4] = 0;
//                    d->qvel[5] = 0;


//                    if (d->time >= 1.5 && d->time <= 1.51) {
//                        d->qvel[1] = 0.7;
//                        d->xfrc_applied[7] = 0.0;
//                        d->xfrc_applied[8] = 0.0;
//                        std::cout << "--- Force Applied ---" << std::endl;
//                    } else {
//                        d->xfrc_applied[7] = 0.0;
//                        d->xfrc_applied[8] = 0.0;
//                    }

                    mj_step(m, d);
                }
                timer.StopTimer();
//                timer.PrintElapsedTime();
            }

//            d->qpos[0] = 0;
//            d->qpos[1] = 0;
//            d->qpos[2] = 0.34;
//            d->qpos[3] = 1;
//            d->qpos[4] = 0;
//            d->qpos[5] = 0;
//            d->qpos[6] = 0;
//
//            d->qvel[0] = 0;
//            d->qvel[1] = 0;
//            d->qvel[2] = 0;
//            d->qvel[3] = 0;
//            d->qvel[4] = 0;
//            d->qvel[5] = 0;

            constexpr double FRAME_RATE = 1.0/60.0;
            const int runs_per_frame = std::floor(FRAME_RATE/dt);
            if (!(run_num_ % runs_per_frame) || pause) {

                // get framebuffer viewport
                mjrRect viewport = {0, 0, 0, 0};
                glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

                // update scene and render
                mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
                UpdateTrajViz();
//            for (int i = 0; i < d->ncon; i++) {
//                mjv_addGeoms(m, d, &opt, NULL, mjCAT_DECOR, &scn);
//                mjtNum force[6];
//                mj_contactForce(m, d, i, force);
//                mjtNum force3[3];
//                force3[0] = .005; //force[0];
//                force3[1] = .005;//force[1];
//                force3[2] = .01*std::sqrt(std::pow(force[0],2) + std::pow(force[1],2) + std::pow(force[2],2));//force[2];
//
//                const float color[4] = {0.5, 0.75, 0.5, 1};
//
//                scn.ngeom += 1;
//                mjv_initGeom(&scn.geoms[scn.ngeom - 1], mjtGeom::mjGEOM_ARROW,
//                             force3, d->contact[i].pos, nullptr, color); //d->contact[i].frame
//            }
                mjr_render(viewport, &scn, &con);

                // swap OpenGL buffers (blocking call due to v-sync)
                glfwSwapBuffers(window);

            }

            // process pending GUI events, call GLFW callbacks
            glfwPollEvents();

            steady_clock::time_point t2 = steady_clock::now();
            if (pause) {
                time_span = duration_cast<duration<double>>(t2 - t2);
            } else {
                time_span = duration_cast<duration<double>>(t2 - t1);
            }
        }

        run_num_++;


        if (glfwWindowShouldClose(window)) {
            //free visualization storage
            mjv_freeScene(&scn);
            mjr_freeContext(&con);

            // free MuJoCo model and data
//            mj_deleteData(d);
//            mj_deleteModel(m);
        }
    }

    void SimpleSimulation::UpdateTrajViz() {
        if (!traj_viz_.empty()) {
//            scn.ngeom = 0;
            const mjtNum mat[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
            const float rgba_va[4] = {1, 1, 0, 1};
            const mjtNum size[3] = {1, 1, 1};
            const mjtNum pos[3] = {1, 1, 1};
            mjtNum from[3];
            mjtNum to[3];
            mjtNum box_pos[3];
            mjtNum box_size[3];

            for (int ee = 0; ee < traj_viz_.size(); ee++) { // end effectors
                for (int node = 0; node < traj_viz_.at(ee).size() - 1; node++) {   // node
                    for (int coord = 0; coord < 3; coord++) {   // coordinate
                        from[coord] = traj_viz_.at(ee).at(node)(coord);
                        to[coord] = traj_viz_.at(ee).at(node+1)(coord);
                    }

                    scn.ngeom += 1;
                    mjv_initGeom(&scn.geoms[scn.ngeom - 1],
                                 mjtGeom::mjGEOM_CAPSULE,
                                 size,
                                 pos,
                                 nullptr,
                                 rgba_va);

                    mjv_connector(&scn.geoms[scn.ngeom - 1],
                                  mjtGeom::mjGEOM_CAPSULE,
                                  0.001,
                                  from,
                                  to);
                }

                if (ee < 4) {
                    for (int i = 0; i < 2; i++) {
                        box_size[i] = box_sides_(i)/2;
                        box_pos[i] = traj_viz_.at(4).at(0)(i) + box_centers_.at(ee)(i);
                    }
                    box_pos[2] = 0.05;
                    box_size[2] = 0.05;

                    const float box_color[4] = {0.5, 0, 0.5, 0.1};

                    scn.ngeom += 1;
                    mjv_initGeom(&scn.geoms[scn.ngeom - 1],
                                 mjtGeom::mjGEOM_BOX,
                                 box_size,
                                 box_pos,
                                 nullptr,
                                 box_color);
                }
            }
        }
    }

    void SimpleSimulation::GetTrajViz(const std::vector<std::vector<Eigen::Vector3d>>& traj_viz,
                                const Eigen::Vector2d& box_sides,
                                const std::vector<Eigen::Vector2d>& box_centers) {
        traj_viz_ = traj_viz;
        box_sides_ = box_sides;
        box_centers_ = box_centers;
    }

    mjtNum* SimpleSimulation::GetControlPointer() {
        return d->ctrl;
    }

    mjData* SimpleSimulation::GetDataPointer() {
        return d;
    }

    const mjModel* SimpleSimulation::GetModelPointer() const {
        return m;
    }

    void SimpleSimulation::PauseSim() {
        pause = true;
    }

} // simulation