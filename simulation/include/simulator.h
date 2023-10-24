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
    extern mjuiState uistate_;

    // mouse interaction
    extern bool button_left_;
    extern bool button_middle_;
    extern bool button_right_;
    extern double lastx_;
    extern double lasty_;

    extern bool run;

    enum {
        // left ui
        SECT_MPC = 0,
        SECT_SIMULATION,
        SECT_PHYSICS,
        SECT_RENDERING,
        NSECT0,

        // right ui
        SECT_JOINT = 0,
        SECT_CONTROL,
        NSECT1
    };

    const double syncMisalign = 0.1;        // maximum mis-alignment before re-sync (simulation seconds)
    const double simRefreshFraction =0.7;


    int UiPredicate(int category, void* userdata);

    namespace util {
        // Utility functions
        template <std::size_t N>
        static inline char* strcpy_arr(char (&dest)[N], const char* src) {
            {
                std::size_t i = 0;
                for (; src[i] && i < N - 1; ++i) {
                    dest[i] = src[i];
                }
                dest[i] = '\0';
            }
            return &dest[0];
        }

        template <std::size_t N>
        static inline int sprintf_arr(char (&dest)[N], const char* format, ...) {
            std::va_list vargs;
            va_start(vargs, format);
            int retval = std::vsnprintf(dest, N, format, vargs);
            va_end(vargs);
            return retval;
        }
    } // util

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
        Simulator(std::unique_ptr<SimulationRobot>& robot, bool viz, double viz_rate);

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

        static constexpr int kMaxFilenameLength = 1000;

        // ----------------------- Member Variables ----------------------- //
        std::recursive_mutex mtx_;
        bool speed_changed;

        int real_time_index = 0;
        int refresh_rate = 60;
        float measured_slowdown = 1.0;

    private:
        // ----------------------- GLFW Interaction Functions ----------------------- //
        static void MouseMove(GLFWwindow *window, double xpos, double ypos);
        static void MouseButton(GLFWwindow* window, int button, int act, int mods);
        static void Scroll(GLFWwindow* window, double xoffset, double yoffset);
        static mjtButton TranslateButton(const int button);
        static void UpdateMJuiState(GLFWwindow* window, mjuiState& state_);

        // -----------------------  ----------------------- //
        void InitSensorFig();
        void UpdateSensor(const mjModel* m, const mjData* d);
        void ShowSensor(mjrRect);

        void UpdateInfoText(const mjModel* model, const mjData* data,
                            char (&title)[kMaxFilenameLength],
                            char (&content)[kMaxFilenameLength]);

        void MakeUISections(const mjModel* m, const mjData* d);
        void MakeControlSection(int oldstate);
        void MakeJointSection(int oldstate);

        void Render();

        void UiModify(mjUI* ui, mjuiState* state, mjrContext* con);
        void UiLayout(mjuiState* state);

        static void UIEvent(mjuiState* state);

        // update an entire section of ui0
        void mjui0_update_section(int section);

        // ----------------------- Member Variables ----------------------- //
        GLFWwindow* window_;
// TODO: Change color
        bool viz_;
        double visualization_rate_;
        bool info = true;
        bool profiler = true;
        bool sensor = false;
        bool update_sensor = false;

        std::unique_ptr<SimulationRobot> robot_;

        // Joint and control info (for viz)
        std::vector<int> jnt_type_;
        std::vector<int> jnt_group_;
        std::vector<int> jnt_qposadr_;
        std::vector<std::optional<std::pair<mjtNum, mjtNum>>> jnt_range_;
        std::vector<std::string> jnt_names_;

        std::vector<int> actuator_group_;
        std::vector<std::optional<std::pair<mjtNum, mjtNum>>> actuator_ctrlrange_;
        std::vector<std::string> actuator_names_;

        std::vector<mjtNum> history_;  // history buffer (nhistory x state_size)

        // Graphics
        mjrContext context_;
        mjvFigure figsensor_ = {};
        double fps_ = 0;

        int spacing = 0;
        int color = 0;

        struct {
            std::optional<std::string> save_xml;
            std::optional<std::string> save_mjb;
            std::optional<std::string> print_model;
            std::optional<std::string> print_data;
            bool reset;
            bool align;
            bool copy_pose;
            bool load_from_history;
            bool load_key;
            bool save_key;
            bool zero_ctrl;
            int newperturb;
            bool select;
            mjuiState select_state;
            bool ui_update_simulation;
            bool ui_update_physics;
            bool ui_update_rendering;
            bool ui_update_joint;
            bool ui_update_ctrl;
        } pending_ = {};

        // info strings
        char info_title[Simulator::kMaxFilenameLength] = {0};
        char info_content[Simulator::kMaxFilenameLength] = {0};

        // simulation section of UI
        const mjuiDef def_simulation[14] = {
                {mjITEM_SECTION,   "Simulation",    1, nullptr,              "AS"},
                {mjITEM_RADIO,     "",              5, &run,           "Pause\nRun"},
                {mjITEM_BUTTON,    "Reset",         2, nullptr,              " #259"},
                {mjITEM_CHECKINT,  "Info",          2, &this->info,       " #291"},
                {mjITEM_BUTTON,    "Reload",        5, nullptr,              "CL"},
                {mjITEM_SEPARATOR, "History",       1},
                {mjITEM_END}
        };

        bool ui0_enable = true;
        bool ui1_enable = true;

        mjUI ui0 = {};
        mjUI ui1 = {};

    }; // Simulator

    void SimulateLoop(Simulator* sim);

}   // simulator


#endif //BILEVEL_GAIT_GEN_SIMULATOR_H
