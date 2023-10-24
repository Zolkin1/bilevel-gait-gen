//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#include <iostream>
#include <functional>
#include <thread>
#include <chrono>
#include <mutex>

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

    bool run = true;

    mjuiState uistate_;



    // logarithmically spaced real-time slow-down coefficients (percent)
    static constexpr float percentRealTime[] = {
            100, 80, 66,  50,  40, 33,  25,  20, 16,  13,
            10,  8,  6.6, 5.0, 4,  3.3, 2.5, 2,  1.6, 1.3,
            1,  .8, .66, .5,  .4, .33, .25, .2, .16, .13,
            .1 };

    Simulator::Simulator(std::unique_ptr<SimulationRobot>& robot, bool viz, double viz_rate) : robot_(std::move(robot)) {
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
        viz_ = viz;
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

        // Setup things for the UI
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window_, &viewport.width, &viewport.height);
        uistate_.nrect = 1;
        uistate_.rect[0].width = viewport.width;
        uistate_.rect[0].height = viewport.height;

        this->ui0.spacing = mjui_themeSpacing(this->spacing);
        this->ui0.color = mjui_themeColor(this->color);
        this->ui0.predicate = &UiPredicate;
        this->ui0.rectid = 1;
        this->ui0.auxid = 0;

        this->ui1.spacing = mjui_themeSpacing(this->spacing);
        this->ui1.color = mjui_themeColor(this->color);
        this->ui1.predicate = UiPredicate;
        this->ui1.rectid = 2;
        this->ui1.auxid = 1;

        mjui_add(&ui0, def_simulation);

        // Set up joint and control info
        jnt_range_.clear();
        jnt_range_.reserve(model_->njnt);
        for (int i = 0; i < model_->njnt; ++i) {
            if (model_->jnt_limited[i]) {
                jnt_range_.push_back(
                        std::make_pair(model_->jnt_range[2 * i], model_->jnt_range[2 * i + 1]));
            } else {
                jnt_range_.push_back(std::nullopt);
            }
        }

        jnt_names_.clear();
        jnt_names_.reserve(model_->njnt);
        for (int i = 0; i < model_->njnt; ++i) {
            jnt_names_.emplace_back(model_->names + model_->name_jntadr[i]);
        }

        actuator_group_.resize(model_->nu);
        std::memcpy(actuator_group_.data(), model_->actuator_group,
                    sizeof(model_->actuator_group[0]) * model_->nu);

        actuator_ctrlrange_.clear();
        actuator_ctrlrange_.reserve(model_->nu);
        for (int i = 0; i < model_->nu; ++i) {
            if (model_->actuator_ctrllimited[i]) {
                actuator_ctrlrange_.push_back(std::make_pair(
                        model_->actuator_ctrlrange[2 * i], model_->actuator_ctrlrange[2 * i + 1]));
            } else {
                actuator_ctrlrange_.push_back(std::nullopt);
            }
        }

        actuator_names_.clear();
        actuator_names_.reserve(model_->nu);
        for (int i = 0; i < model_->nu; ++i) {
            actuator_names_.emplace_back(model_->names + model_->name_actuatoradr[i]);
        }

        // Make UI sections
        MakeUISections(model_, data_);

        // Setup userdata
        uistate_.userdata = this;

        // Set up the UI
        InitSensorFig();

        //TODO: Set the control callback function (not the mujoco one)

    }

    void Simulator::RunSimulator() {

        // Create a simulation thread
        std::thread simulationThreadHandle(&SimulateLoop, this);

        while (!glfwWindowShouldClose(window_)) {
            // advance interactive simulation for 1/60 sec
            //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
            //  this loop will finish on time for the next frame to be rendered at 60 fps.
            //  Otherwise add a cpu timer and exit this loop when it is time to render.
            mjtNum simstart = data_->time;

            // take the mutex and prepare for the render
            std::unique_lock<std::recursive_mutex> lock(mtx_);

            // update the visual scene
            mjv_updateScene(model_, data_, &opt_, nullptr, &cam_, mjCAT_ALL, &scene_);

            // update sensor
            UpdateSensor(model_, data_);

            // release the mutex
            lock.unlock();

            Render();

            // Render while simulating
//            mjrRect rect = this->uistate_.rect[3];
//            mjrRect smallrect = rect;

//            ShowSensor(viewport);


            // get framebuffer viewport
            //mjrRect viewport = {0, 0, 0, 0};
//            glfwGetFramebufferSize(window_, &viewport.width, &viewport.height);
//
//            // update scene and render
//            //mjv_updateScene(model_, data_, &opt_, NULL, &cam_, mjCAT_ALL, &scene_);
//            //mjr_render(viewport, &scene_, &context_);
//            mjr_render(viewport, &scene_, &context_);
//
//            // swap OpenGL buffers (blocking call due to v-sync)
//            glfwSwapBuffers(window_);

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

    // ----------------------- UI Functions ----------------------- //
    void Simulator::Render() {
        // Allow for window re-sizing
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window_, &viewport.width, &viewport.height);
//        mjr_resizeOffscreen(viewport.width, viewport.height, &context_);
        uistate_.rect[0].width = viewport.width;
        uistate_.rect[0].height = viewport.height;


        // update rendering context buffer size if required
        UiModify(&this->ui0, &uistate_, &context_);
        UiModify(&this->ui1, &uistate_, &context_);

        // get 3D rectangle and reduced for profiler
        mjrRect rect = uistate_.rect[3];
        mjrRect smallrect = rect;
        smallrect.width = rect.width - rect.width/4;

        // update UI sections from last sync
        if (pending_.ui_update_rendering) {
            if (this->ui0_enable && this->ui0.sect[SECT_MPC].state) {
                mjui0_update_section(SECT_MPC);
            }
            pending_.ui_update_rendering = false;
        }

        if (pending_.ui_update_simulation) {
            if (this->ui0_enable && this->ui0.sect[SECT_SIMULATION].state) {
                mjui0_update_section(SECT_SIMULATION);
            }
            pending_.ui_update_simulation = false;
        }

        if (pending_.ui_update_physics) {
            if (this->ui0_enable && this->ui0.sect[SECT_PHYSICS].state) {
                mjui0_update_section(SECT_PHYSICS);
            }
            pending_.ui_update_physics = false;
        }

        if (pending_.ui_update_rendering) {
            if (this->ui0_enable && this->ui0.sect[SECT_RENDERING].state) {
                mjui0_update_section(SECT_RENDERING);
            }
            pending_.ui_update_rendering = false;
        }

        if (pending_.ui_update_joint) {
            if (this->ui1_enable && this->ui1.sect[SECT_JOINT].state) {
                mjui_update(SECT_JOINT, -1, &this->ui1, &uistate_, &context_);
            }
            pending_.ui_update_joint = false;
        }

        if (pending_.ui_update_ctrl) {
            if (this->ui1_enable && this->ui1.sect[SECT_CONTROL].state) {
                mjui_update(SECT_CONTROL, -1, &this->ui1, &uistate_, &context_);
            }
            pending_.ui_update_ctrl = false;
        }

        // Sync
        // Get joint info
        for (int i = 0; i < model_->njnt; ++i) {
            std::optional<std::pair<mjtNum, mjtNum>> range;
            if (model_->jnt_limited[i]) {
                range.emplace(model_->jnt_range[2*i], model_->jnt_range[2*i + 1]);
            }
            if (jnt_range_[i] != range) {
                pending_.ui_update_joint = true;
                jnt_range_[i].swap(range);
            }
        }

        // Get control info
        for (int i = 0; i < model_->nu; ++i) {
            std::optional<std::pair<mjtNum, mjtNum>> range;
            if (model_->actuator_ctrllimited[i]) {
                range.emplace(model_->actuator_ctrlrange[2*i], model_->actuator_ctrlrange[2*i + 1]);
            }
            if (actuator_ctrlrange_[i] != range) {
                pending_.ui_update_ctrl = true;
                actuator_ctrlrange_[i].swap(range);
            }
        }

        // Check the reset
        if (pending_.reset) {
            mj_resetData(model_, data_);
            mj_forward(model_, data_);
            update_sensor = true;
            pending_.ui_update_simulation = true;
            pending_.reset = false;
        }

        // update
        mjv_updateScene(model_, data_, &opt_, nullptr, &cam_, mjCAT_ALL, &scene_);

        if (this->info) {
            UpdateInfoText(model_, data_, this->info_title, this->info_content);
        }

        // render scene
        mjr_render(rect, &scene_, &context_);

        // get desired and actual percent-of-real-time
        float desiredRealtime = percentRealTime[this->real_time_index];
        float actualRealtime = 100 / this->measured_slowdown;

        // if running, check for misalignment of more than 10%
        float realtime_offset = mju_abs(actualRealtime - desiredRealtime);
        bool misaligned = realtime_offset > 0.1 * desiredRealtime;

        // make realtime overlay label
        char rtlabel[30] = {'\0'};
        if (desiredRealtime != 100.0 || misaligned) {
            // print desired realtime
            int labelsize = std::snprintf(rtlabel, sizeof(rtlabel), "%g%%", desiredRealtime);

            // if misaligned, append to label
            if (misaligned) {
                std::snprintf(rtlabel+labelsize, sizeof(rtlabel)-labelsize, " (%-4.1f%%)", actualRealtime);
            }
        }

        // show real-time overlay
        if (rtlabel[0]) {
            mjr_overlay(mjFONT_BIG, mjGRID_TOPLEFT, smallrect, rtlabel, nullptr, &context_);
        }

        // show ui 0
        if (this->ui0_enable) {
            mjui_render(&this->ui0, &uistate_, &context_);
        }

        // show ui 1
        if (this->ui1_enable) {
            mjui_render(&this->ui1, &uistate_, &context_);
        }

        // show info
        if (this->info) {
            mjr_overlay(mjFONT_NORMAL, mjGRID_BOTTOMLEFT, rect, info_title, info_content,
                        &context_);
        }

        // show profiler
        if (this->profiler) {
            //ShowProfiler(this, rect);
        }

        // show sensor
        if (this->sensor) {
            ShowSensor(smallrect);
        }

        // finalize
        glfwSwapBuffers(window_);
        //this->platform_ui->SwapBuffers();
    }

    void Simulator::UiModify(mjUI* ui, mjuiState* state, mjrContext* con) {
        mjui_resize(ui, con);
        mjr_addAux(ui->auxid, ui->width, ui->maxheight, ui->spacing.samples, con);
        UiLayout(state);
        mjui_update(-1, -1, ui, state, con);
    }

    // set window layout
    void Simulator::UiLayout(mjuiState* state) {
        mjrRect* rect = state->rect;

        // set number of rectangles
        state->nrect = 4;

        // rect 1: UI 0
        rect[1].left = 0;
        rect[1].width = ui0_enable ? ui0.width : 0;
        rect[1].bottom = 0;
        rect[1].height = rect[0].height;

        // rect 2: UI 1
        rect[2].width = ui1_enable ? ui1.width : 0;
        rect[2].left = mjMAX(0, rect[0].width - rect[2].width);
        rect[2].bottom = 0;
        rect[2].height = rect[0].height;

        // rect 3: 3D plot (everything else is an overlay)
        rect[3].left = rect[1].width;
        rect[3].width = mjMAX(0, rect[0].width - rect[1].width - rect[2].width);
        rect[3].bottom = 0;
        rect[3].height = rect[0].height;
    }

    // determine enable/disable item state given category
    int UiPredicate(int category, void* userdata) {

        switch (category) {
            case 2:                 // require model
                return 1;

            case 3:                 // require model and nkey
                return 0;

            case 4:                 // require model and paused
                return !run;

            case 5:                 // require model and fully managed mode
                return 0;

            default:
                return 1;
        }
    }

    void Simulator::mjui0_update_section(int section) {
        mjui_update(section, -1, &ui0, &uistate_, &context_);
    }

    void Simulator::InitSensorFig() {
        mjv_defaultFigure(&figsensor_);

        figsensor_.figurergba[3] = 0.5f;

        // set flags
        figsensor_.flg_extend = 1;
        figsensor_.flg_barplot = 1;
        figsensor_.flg_symmetric = 1;

        // title
        util::strcpy_arr(figsensor_.title, "Sensor data");

        // y-tick number format
        util::strcpy_arr(figsensor_.yformat, "%.0f");

        // grid size
        figsensor_.gridsize[0] = 2;
        figsensor_.gridsize[1] = 3;

        // minimum range
        figsensor_.range[0][0] = 0;
        figsensor_.range[0][1] = 0;
        figsensor_.range[1][0] = -1;
        figsensor_.range[1][1] = 1;
    }

    void Simulator::UpdateSensor(const mjModel* m, const mjData* d) {
        static const int maxline = 10;

        // clear linepnt
        for (int i=0; i<maxline; i++) {
            figsensor_.linepnt[i] = 0;
        }

        // start with line 0
        int lineid = 0;

        // loop over sensors
        for (int n=0; n<m->nsensor; n++) {
            // go to next line if type is different
            if (n>0 && m->sensor_type[n]!=m->sensor_type[n-1]) {
                lineid = mjMIN(lineid+1, maxline-1);
            }

            // get info about this sensor
            mjtNum cutoff = (m->sensor_cutoff[n]>0 ? m->sensor_cutoff[n] : 1);
            int adr = m->sensor_adr[n];
            int dim = m->sensor_dim[n];

            // data pointer in line
            int p = figsensor_.linepnt[lineid];

            // fill in data for this sensor
            for (int i=0; i<dim; i++) {
                // check size
                if ((p+2*i)>=mjMAXLINEPNT/2) {
                    break;
                }

                // x
                figsensor_.linedata[lineid][2*p+4*i] = adr+i;
                figsensor_.linedata[lineid][2*p+4*i+2] = adr+i;

                // y
                figsensor_.linedata[lineid][2*p+4*i+1] = 0;
                figsensor_.linedata[lineid][2*p+4*i+3] = d->sensordata[adr+i]/cutoff;
            }

            // update linepnt
            figsensor_.linepnt[lineid] = mjMIN(mjMAXLINEPNT-1, figsensor_.linepnt[lineid]+2*dim);
        }
    }

    void Simulator::UpdateInfoText(const mjModel* model, const mjData* data,
                        char (&title)[kMaxFilenameLength],
                        char (&content)[kMaxFilenameLength]) {
        char tmp[20];

        // number of islands with statistics
        int nisland = mjMIN(data_->solver_nisland, mjNISLAND);

        // compute solver error (maximum over islands)
        mjtNum solerr = 0;
        for (int i=0; i < nisland; i++) {
            mjtNum solerr_i = 0;
            if (data_->solver_niter[i]) {
                int ind = mjMIN(data_->solver_niter[i], mjNSOLVER) - 1;
                const mjSolverStat* stat = data_->solver + i*mjNSOLVER + ind;
                solerr_i = mju_min(stat->improvement, stat->gradient);
                if (solerr_i==0) {
                    solerr_i = mju_max(stat->improvement, stat->gradient);
                }
            }
            solerr = mju_max(solerr, solerr_i);
        }
        solerr = mju_log10(mju_max(mjMINVAL, solerr));

        // format FPS text
        char fps[10];
        if (fps_ < 1) {
            util::sprintf_arr(fps, "%0.1f ", fps_);
        } else {
            util::sprintf_arr(fps, "%.0f ", fps_);
        }

        // total iterations of all islands with statistics
        int solver_niter = 0;
        for (int i=0; i < nisland; i++) {
            solver_niter += data_->solver_niter[i];
        }

        // prepare info text
        util::strcpy_arr(title, "Time\nSize\nCPU\nSolver   \nFPS\nMemory");
        util::sprintf_arr(content,
                         "%-9.3f\n%d  (%d con)\n%.3f\n%.1f  (%d it)\n%s\n%.2g of %s",
                         data_->time,
                         data_->nefc, data_->ncon,
                         run ?
                         data_->timer[mjTIMER_STEP].duration / mjMAX(1, data_->timer[mjTIMER_STEP].number) :
                         data_->timer[mjTIMER_FORWARD].duration / mjMAX(1, data_->timer[mjTIMER_FORWARD].number),
                         solerr, solver_niter,
                         fps,
                         data_->maxuse_arena/(double)(data_->narena),
                         mju_writeNumBytes(data_->narena));
    }

    // show sensor figure
    void Simulator::ShowSensor(mjrRect rect) {
        // constant width with and without profiler
        int width = rect.width/4; //sim->profiler ? rect.width/3 : rect.width/4;

        // render figure on the right
        mjrRect viewport = {
                rect.left + rect.width - width,
                rect.bottom,
                width,
                rect.height/3
        };
        mjr_figure(viewport, &figsensor_, &context_);
    }

    void Simulator::MakeUISections(const mjModel* m, const mjData* d) {
        // get section open-close state, UI 0
        int oldstate0[NSECT0];
        for (int i=0; i<NSECT0; i++) {
            oldstate0[i] = 0;
            if (ui0.nsect>i) {
                oldstate0[i] = ui0.sect[i].state;
            }
        }

        // get section open-close state, UI 1
        int oldstate1[NSECT1];
        for (int i=0; i<NSECT1; i++) {
            oldstate1[i] = 0;
            if (ui1.nsect>i) {
                oldstate1[i] = ui1.sect[i].state;
            }
        }

        // clear model-dependent sections of UI
        ui0.nsect = SECT_PHYSICS;
        ui1.nsect = 0;

        // make
        //MakeRenderingSection(sim, m, oldstate0[SECT_RENDERING]);
        //MakeVisualizationSection(sim, m, oldstate0[SECT_VISUALIZATION]);
        MakeJointSection(oldstate1[SECT_JOINT]);
        MakeControlSection(oldstate1[SECT_CONTROL]);
    }

    void Simulator::MakeJointSection(int oldstate) {
        mjuiDef defJoint[] = {
                {mjITEM_SECTION, "Joint", oldstate, nullptr, "AJ"},
                {mjITEM_END}
        };
        mjuiDef defSlider[] = {
                {mjITEM_SLIDERNUM, "", 2, nullptr, "0 1"},
                {mjITEM_END}
        };

        // add section
        mjui_add(&ui1, defJoint);
        defSlider[0].state = 4;

        // add scalar joints, exit if UI limit reached
        int itemcnt = 0;
        for (int i=0; i < jnt_type_.size() && itemcnt<mjMAXUIITEM; i++)
            if ((jnt_type_[i]==mjJNT_HINGE || jnt_type_[i]==mjJNT_SLIDE)) {
                // skip if joint group is disabled
                if (!opt_.jointgroup[mjMAX(0, mjMIN(mjNGROUP-1, jnt_group_[i]))]) {
                    continue;
                }

                // set data and name
                defSlider[0].pdata = &data_->qpos[model_->jnt_qposadr[i]];

                if (!jnt_names_[i].empty()) {
                    util::strcpy_arr(defSlider[0].name, jnt_names_[i].c_str());
                } else {
                    util::sprintf_arr(defSlider[0].name, "joint %d", i);
                }

                // set range
                if (jnt_range_[i].has_value())
                    util::sprintf_arr(defSlider[0].other, "%.4g %.4g",
                                     jnt_range_[i]->first, jnt_range_[i]->second);
                else if (jnt_type_[i]==mjJNT_SLIDE) {
                    util::strcpy_arr(defSlider[0].other, "-1 1");
                } else {
                    util::strcpy_arr(defSlider[0].other, "-3.1416 3.1416");
                }

                // add and count
                mjui_add(&ui1, defSlider);
                itemcnt++;
            }
    }

    // make control section of UI
    void Simulator::MakeControlSection(int oldstate) {
        mjuiDef defControl[] = {
                {mjITEM_SECTION, "Control", oldstate, nullptr, "AC"},
                {mjITEM_BUTTON,  "Clear all", 2},
                {mjITEM_END}
        };
        mjuiDef defSlider[] = {
                {mjITEM_SLIDERNUM, "", 2, nullptr, "0 1"},
                {mjITEM_END}
        };

        // add section
        mjui_add(&ui1, defControl);
        defSlider[0].state = 2;

        // add controls, exit if UI limit reached (Clear button already added)
        int itemcnt = 1;
        for (int i=0; i < actuator_ctrlrange_.size() && itemcnt<mjMAXUIITEM; i++) {
            // skip if actuator group is disabled
            if (!opt_.actuatorgroup[mjMAX(0, mjMIN(mjNGROUP-1, actuator_group_[i]))]) {
                continue;
            }

            // set data and name
            defSlider[0].pdata = &data_->ctrl[i];

            if (!actuator_names_[i].empty()) {
                util::strcpy_arr(defSlider[0].name, actuator_names_[i].c_str());
            } else {
                util::sprintf_arr(defSlider[0].name, "control %d", i);
            }

            // set range
            if (actuator_ctrlrange_[i].has_value())
                util::sprintf_arr(defSlider[0].other, "%.4g %.4g",
                                 actuator_ctrlrange_[i]->first, actuator_ctrlrange_[i]->second);
            else {
                util::strcpy_arr(defSlider[0].other, "-1 1");
            }

            // add and count
            mjui_add(&ui1, defSlider);
            itemcnt++;
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

    void Simulator::UpdateMJuiState(GLFWwindow* window, mjuiState& state_) {
        // mouse buttons
        state_.right = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
        state_.left = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        state_.middle = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;

        // keyboard modifiers
//        state_.control = IsCtrlKeyPressed();
//        state_.shift = IsShiftKeyPressed();
//        state_.alt = IsAltKeyPressed();

        // swap left and right if Alt
//        if (state_.alt) {
//            int tmp = state_.left;
//            state_.left = state_.right;
//            state_.right = tmp;
//        }

        // get mouse position, scale by buffer-to-window ratio'
        double x, y;
        glfwGetCursorPos(window, &lastx_, &lasty_);

        int w1, h1, w2, h2;
        glfwGetFramebufferSize(window, &w1, &h1);
        glfwGetWindowSize(window, &w2, &h2);
        const double buffer_window_ratio = static_cast<double>(w1) / w2;
        x *= buffer_window_ratio;
        y *= buffer_window_ratio;

        // invert y to match OpenGL convention
        y = state_.rect[0].height - y;

        // save
        state_.dx = x - state_.x;
        state_.dy = y - state_.y;
        state_.x = x;
        state_.y = y;

        // find mouse rectangle
        state_.mouserect = mjr_findRect(mju_round(x), mju_round(y), state_.nrect - 1, state_.rect + 1) + 1;
    }

    // mouse button callback
    void Simulator::MouseButton(GLFWwindow* window, int button, int act, int mods) {
//        // update button state
//        mjtButton mj_button = TranslateButton(button);
//        uistate_.right = 1;
//        // update mouse position
//        glfwGetCursorPos(window, &lastx_, &lasty_);
//        int w1, h1, w2, h2;
//        glfwGetFramebufferSize(window, &w1, &h1);
//        glfwGetWindowSize(window, &w2, &h2);
//        const double buffer_window_ratio = static_cast<double>(w1)/w2;
//        //lastx_ *= buffer_window_ratio;
//        //lasty_ *= buffer_window_ratio;
//
//        // find mouse rectangle
//        uistate_.mouserect = mjr_findRect(mju_round(lastx_), mju_round(lasty_), uistate_.nrect-1, uistate_.rect+1) + 1;
//        uistate_.button = mj_button;
//
//        double now = std::chrono::duration<double>(
//                std::chrono::steady_clock::now().time_since_epoch()).count();
//        uistate_.buttontime = now;
//        uistate_.type = mjEVENT_RELEASE;
//
//        uistate_.x = lastx_;
//        uistate_.y = lasty_;
//
//        UIEvent(&uistate_);

        mjtButton mj_button = TranslateButton(button);

        // update state
        UpdateMJuiState(window, uistate_);

        // swap left and right if Alt
        if (uistate_.alt) {
            if (mj_button == mjBUTTON_LEFT) {
                mj_button = mjBUTTON_RIGHT;
            } else if (mj_button == mjBUTTON_RIGHT) {
                mj_button = mjBUTTON_LEFT;
            }
        }

        // press
        if (GLFW_PRESS) {
            double now = std::chrono::duration<double>(
                    std::chrono::steady_clock::now().time_since_epoch()).count();

            // detect doubleclick: 250 ms
            if (mj_button == uistate_.button && now - uistate_.buttontime < 0.25) {
                uistate_.doubleclick = 1;
            } else {
                uistate_.doubleclick = 0;
            }

            // set info
            uistate_.type = mjEVENT_PRESS;
            uistate_.button = mj_button;
            uistate_.buttontime = now;

            // start dragging
            if (uistate_.mouserect) {
                uistate_.dragbutton = uistate_.button;
                uistate_.dragrect = uistate_.mouserect;
            }
        }

            // release
        else {
            uistate_.type = mjEVENT_RELEASE;
        }

        // application-specific processing
        UIEvent(&uistate_);

        // stop dragging after application processing
        if (uistate_.type == mjEVENT_RELEASE) {
            uistate_.dragrect = 0;
            uistate_.dragbutton = 0;
        }
    }

    mjtButton Simulator::TranslateButton(const int button) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            return mjBUTTON_LEFT;
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            return mjBUTTON_RIGHT;
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
            return mjBUTTON_MIDDLE;
        }
        return mjBUTTON_NONE;
    }

    // scroll callback
    void Simulator::Scroll(GLFWwindow* window, double xoffset, double yoffset) {
        // emulate vertical mouse motion = 5% of window height
        mjv_moveCamera(model_, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scene_, &cam_);
    }

    // ----------------------- Simulate Loop ----------------------- //
    void SimulateLoop(Simulator* sim) {

        // cpu-sim syncronization point
        std::chrono::time_point<std::chrono::steady_clock> syncCPU;
        mjtNum syncSim = 0;

        // run always
        while (1) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            {
                // lock the sim mutex
                const std::unique_lock<std::recursive_mutex> lock(sim->mtx_);

                bool stepped = false;

                // record cpu time at start of iteration
                const auto startCPU = std::chrono::steady_clock::now();

                // elapsed CPU and simulation time since last sync
                const auto elapsedCPU = startCPU - syncCPU;
                double elapsedSim = data_->time - syncSim;

                // requested slow-down factor
                double slowdown = 100 / percentRealTime[sim->real_time_index];

                // misalignment condition: distance from target sim time is bigger than syncmisalign
                bool misaligned =
                        mju_abs(std::chrono::duration<double>(elapsedCPU).count()/slowdown - elapsedSim) > syncMisalign;

                // out-of-sync (for any reason): reset sync times, step
                if (elapsedSim < 0 || elapsedCPU.count() < 0 || syncCPU.time_since_epoch().count() == 0 ||
                    misaligned || sim->speed_changed) {
                    // re-sync
                    syncCPU = startCPU;
                    syncSim = data_->time;
                    sim->speed_changed = false;

                    // run single step, let next iteration deal with timing
                    mj_step(model_, data_);
                    stepped = true;
                }

                // in-sync: step until ahead of cpu
                else {
                    bool measured = false;
                    mjtNum prevSim = data_->time;

                    double refreshTime = simRefreshFraction/sim->refresh_rate;

                    // step while sim lags behind cpu and within refreshTime
                    while (std::chrono::duration<double>((data_->time - syncSim)*slowdown) < std::chrono::steady_clock::now() - syncCPU &&
                            std::chrono::steady_clock::now() - startCPU < std::chrono::duration<double>(refreshTime)) {
                        // measure slowdown before first step
                        if (!measured && elapsedSim) {
                            sim->measured_slowdown =
                                    std::chrono::duration<double>(elapsedCPU).count() / elapsedSim;
                            measured = true;
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
//                if (stepped) {
//                    sim.AddToHistory();
//                }


            }  // release std::lock_guard<std::mutex>
        }
    }

    // handle UI event
    void Simulator::UIEvent(mjuiState* state) {
        Simulator* sim = static_cast<Simulator*>(state->userdata);

        // call UI 0 if event is directed to it
        if ((state->dragrect==sim->ui0.rectid) ||
            (state->dragrect==0 && state->mouserect==sim->ui0.rectid) ||
            state->type==mjEVENT_KEY) {
                // process UI event
                mjuiItem* it = mjui_event(&sim->ui0, state, &sim->context_);

        // simulation section
        if (it && it->sectionid==SECT_SIMULATION) {
            switch (it->itemid) {
                case 1:             // Reset
                    sim->pending_.reset = true;
                    break;

                case 2:             // Info
                    sim->pending_.copy_pose = true;
                    break;

                case 3:             // Reload
                    //sim->uiloadrequest.fetch_add(1);
                    break;

                case 4:            // History scrubber
                    run = 0;
                    sim->pending_.load_from_history = true;
                    sim->mjui0_update_section(SECT_SIMULATION);
                    break;
            }
        }
        // physics section
//        else if (it && it->sectionid==SECT_PHYSICS && sim->m_) {
//            mjOption* opt = sim->is_passive_ ? &sim->scnstate_.model.opt : &sim->m_->opt;
//
//            // update disable flags in mjOption
//            opt->disableflags = 0;
//            for (int i=0; i<mjNDISABLE; i++)
//                if (sim->disable[i]) {
//                    opt->disableflags |= (1<<i);
//                }
//
//            // update enable flags in mjOption
//            opt->enableflags = 0;
//            for (int i=0; i<mjNENABLE; i++)
//                if (sim->enable[i]) {
//                    opt->enableflags |= (1<<i);
//                }
//            }

            // rendering section
//            else if (it && it->sectionid==SECT_RENDERING) {
//                // set camera in mjvCamera
//                if (sim->camera==0) {
//                    sim->cam.type = mjCAMERA_FREE;
//                } else if (sim->camera==1) {
//                    if (sim->pert.select>0) {
//                        sim->cam.type = mjCAMERA_TRACKING;
//                        sim->cam.trackbodyid = sim->pert.select;
//                        sim->cam.fixedcamid = -1;
//                    } else {
//                        sim->cam.type = mjCAMERA_FREE;
//                        sim->camera = 0;
//                        mjui0_update_section(sim, SECT_RENDERING);
//                    }
//                } else {
//                    sim->cam.type = mjCAMERA_FIXED;
//                    sim->cam.fixedcamid = sim->camera - 2;
//                }
//                // copy camera spec to clipboard (as MJCF element)
//                if (it->itemid == 3) {
//                    CopyCamera(sim);
//                }
//            }

                // visualization section
//            else if (it && it->sectionid==SECT_VISUALIZATION) {
//                if (!mju::strcmp_arr(it->name, "Align")) {
//                    sim->pending_.align = true;
//                }
//            }

            // stop if UI processed event
            if (it!=nullptr || (state->type==mjEVENT_KEY && state->key==0)) {
                return;
            }
        }

        // call UI 1 if event is directed to it
        if ((state->dragrect==sim->ui1.rectid) ||
            (state->dragrect==0 && state->mouserect==sim->ui1.rectid) ||
            state->type==mjEVENT_KEY) {
            // process UI event
            mjuiItem* it = mjui_event(&sim->ui1, state, &sim->context_);

            // control section
            if (it && it->sectionid==SECT_CONTROL) {
                // clear controls
                if (it->itemid==0) {
                    sim->pending_.zero_ctrl = true;
                }
            }

            // stop if UI processed event
            if (it!=nullptr || (state->type==mjEVENT_KEY && state->key==0)) {
                return;
            }
        }

        // shortcut not handled by UI
        if (state->type==mjEVENT_KEY && state->key!=0) {
            switch (state->key) {
                case ' ':                   // Mode
                    if (model_) {
                        run = 1 - run;
                        //sim->pert.active = 0;

                        sim->mjui0_update_section(-1);
                    }
                    break;
            }

            return;
        }

        // 3D scroll
//        if (state->type==mjEVENT_SCROLL && state->mouserect==3) {
//            // emulate vertical mouse motion = 2% of window height
//            if (sim->m_ && !sim->is_passive_) {
//                mjv_moveCamera(sim->m_, mjMOUSE_ZOOM, 0, -zoom_increment*state->sy, &sim->scn, &sim->cam);
//            } else {
//                mjv_moveCameraFromState(
//                        &sim->scnstate_, mjMOUSE_ZOOM, 0, -zoom_increment*state->sy, &sim->scn, &sim->cam);
//            }
//            return;
//        }

        // 3D press
//        if (state->type==mjEVENT_PRESS && state->mouserect==3) {
//            // set perturbation
//            int newperturb = 0;
//            if (state->control && sim->pert.select>0 && (sim->m_ || sim->is_passive_)) {
//                // right: translate;  left: rotate
//                if (state->right) {
//                    newperturb = mjPERT_TRANSLATE;
//                } else if (state->left) {
//                    newperturb = mjPERT_ROTATE;
//                }
//                if (newperturb && !sim->pert.active) {
//                    sim->pending_.newperturb = newperturb;
//                }
//            }
//
//            // handle double-click
//            if (state->doubleclick && (sim->m_ || sim->is_passive_)) {
//                sim->pending_.select = true;
//                std::memcpy(&sim->pending_.select_state, state, sizeof(sim->pending_.select_state));
//
//                // stop perturbation on select
//                sim->pert.active = 0;
//                sim->pending_.newperturb = 0;
//            }
//
//            return;
//        }

        // 3D release
//        if (state->type==mjEVENT_RELEASE && state->dragrect==3 && (sim->m_ || sim->is_passive_)) {
//            // stop perturbation
//            sim->pert.active = 0;
//            sim->pending_.newperturb = 0;
//            return;
//        }

        // 3D move
//        if (state->type==mjEVENT_MOVE && state->dragrect==3 && (sim->m_ || sim->is_passive_)) {
//            // determine action based on mouse button
//            mjtMouse action;
//            if (state->right) {
//                action = state->shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
//            } else if (state->left) {
//                action = state->shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
//            } else {
//                action = mjMOUSE_ZOOM;
//            }
//
//            // move perturb or camera
//            mjrRect r = state->rect[3];
//            if (sim->pert.active) {
//                if (!sim->is_passive_) {
//                    mjv_movePerturb(
//                            sim->m_, sim->d_, action, state->dx / r.height, -state->dy / r.height,
//                            &sim->scn, &sim->pert);
//                } else {
//                    mjv_movePerturbFromState(
//                            &sim->scnstate_, action, state->dx / r.height, -state->dy / r.height,
//                            &sim->scn, &sim->pert);
//                }
//            } else {
//                if (!sim->is_passive_) {
//                    mjv_moveCamera(
//                            sim->m_, action, state->dx / r.height, -state->dy / r.height,
//                            &sim->scn, &sim->cam);
//                } else {
//                    mjv_moveCameraFromState(
//                            &sim->scnstate_, action, state->dx / r.height, -state->dy / r.height,
//                            &sim->scn, &sim->cam);
//                }
//            }
//            return;
//        }

        // Redraw
        if (state->type == mjEVENT_REDRAW) {
            sim->Render();
            return;
        }
    }
} // simulator