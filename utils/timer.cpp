//
// Created by zolkin on 12/29/23.
//

#include "timer.h"

#include <utility>

namespace utils {
    Timer::Timer(std::string  name) : name_(std::move(name)) {
        begin_ = std::chrono::steady_clock::now();
        end_ = std::chrono::steady_clock::now();
        elapsed_time_ = 0;
    }

    void Timer::StartTimer() {
        begin_ = std::chrono::steady_clock::now();
    }

    void Timer::StopTimer() {
        end_ = std::chrono::steady_clock::now();
        elapsed_time_ = std::chrono::duration_cast<std::chrono::milliseconds>(end_ - begin_).count();
        if (elapsed_time_ < 0) {
            std::cerr << "Bad timer. Negative elapsed time." << std::endl;
        }
    }

    long Timer::GetElapsedTimeMilliseconds() const {
        return elapsed_time_;
    }

    double Timer::GetElapsedTimeSeconds() const {
        return static_cast<double>(elapsed_time_)/1000.0;
    }

    void Timer::PrintElapsedTime() const {
        std::cout << "[" << name_ << "] " << "took " << elapsed_time_ << " ms." << std::endl;
    }
}