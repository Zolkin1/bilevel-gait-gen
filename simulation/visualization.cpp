//
// Created by zolkin on 12/12/23.
//

#include <iostream>
#include <chrono>

#include "visualization.h"
#include "simulation_robot.h"

namespace simulation {
    // mouse interaction
    bool button_left = false;
    bool button_middle = false;
    bool button_right =  false;
    double lastx = 0;
    double lasty = 0;

    bool pause = false;
    int step = 0;

    // MuJoCo data structures
    mjModel* m = NULL;                  // MuJoCo model
    mjData* d = NULL;                   // MuJoCo data
    mjvCamera cam;                      // abstract camera
    mjvOption opt;                      // visualization options
    mjvScene scn;                       // abstract scene
    mjrContext con;                     // custom GPU context

    GLFWwindow* window;
// keyboard callback
    void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods) {
        // backspace: reset simulation
        if (act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE) {
            mj_resetData(m, d);
            mj_forward(m, d);
        }

        if (act == GLFW_PRESS && key == GLFW_KEY_SPACE) {
            pause = !pause;
        }

        if (act == GLFW_PRESS && key == GLFW_KEY_S) {
            std::cout << "Step #" << step << ": " << std::endl;
            std::cout << "Configuration: ";
            for (int i = 0; i < m->nq; i++) {
                std::cout << std::setw(15) << std::setprecision(6) << d->qpos[i];
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }


// mouse button callback
    void mouse_button(GLFWwindow* window, int button, int act, int mods) {
        // update button state
        button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
        button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
        button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

        // update mouse position
        glfwGetCursorPos(window, &lastx, &lasty);
    }


// mouse move callback
    void mouse_move(GLFWwindow* window, double xpos, double ypos) {
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
        bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                          glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

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
        mjv_moveCamera(m, action, dx/height, dy/height, &scn, &cam);
    }


// scroll callback
    void scroll(GLFWwindow* window, double xoffset, double yoffset) {
        // emulate vertical mouse motion = 5% of window height
        mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scn, &cam);
    }

    Visualizer::Visualizer(const std::string& robot_xml) {
        // load and compile model
        char error[1000] = "Could not load binary model";
        m = mj_loadXML(robot_xml.c_str(), nullptr, error, 1000);

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
        window = glfwCreateWindow(1200, 900, "Visualizer", NULL, NULL);
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
    }

    void Visualizer::UpdateState(const vector_t& config) {
        assert(config.size() == m->nq);
        mj_step(m,d);
        for (int i = 0; i < m->nq; i++) {
            d->qpos[i] = config(i);
        }

        for (int i = 0; i < m->nv; i++) {
            d->qvel[i] = 0;
            d->qacc[i] = 0;
        }
        step++;
    }

    void Visualizer::UpdateViz(double dt) {
        using namespace std::chrono;
        steady_clock::time_point t1 = steady_clock::now();

        duration<double> time_span(0);
        while(!glfwWindowShouldClose(window) && time_span.count() < dt) {
            // get framebuffer viewport
            mjrRect viewport = {0, 0, 0, 0};
            glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

            // update scene and render
            mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
            UpdateTrajViz();
            mjr_render(viewport, &scn, &con);

            // swap OpenGL buffers (blocking call due to v-sync)
            glfwSwapBuffers(window);

            // process pending GUI events, call GLFW callbacks
            glfwPollEvents();

            steady_clock::time_point t2 = steady_clock::now();
            if (pause) {
                time_span = duration_cast<duration<double>>(t2 - t2);
            } else {
                time_span = duration_cast<duration<double>>(t2 - t1);
            }
        }

        if (glfwWindowShouldClose(window)) {
            //free visualization storage
            mjv_freeScene(&scn);
            mjr_freeContext(&con);

            // free MuJoCo model and data
//            mj_deleteData(d);
//            mj_deleteModel(m);
        }
    }

    const mjModel* Visualizer::GetModel() const {
        return m;
    }

    void Visualizer::UpdateTrajViz() {
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

    void Visualizer::GetTrajViz(const std::vector<std::vector<Eigen::Vector3d>>& traj_viz,
                                const Eigen::Vector2d& box_sides,
                                const std::vector<Eigen::Vector2d>& box_centers) {
        traj_viz_ = traj_viz;
        box_sides_ = box_sides;
        box_centers_ = box_centers;
    }

}