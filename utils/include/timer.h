//
// Created by zolkin on 12/29/23.
//

#ifndef BILEVEL_GAIT_GEN_TIMER_H
#define BILEVEL_GAIT_GEN_TIMER_H

#include <chrono>
#include <string>
#include <iostream>

namespace utils {
    class Timer {
    public:
        Timer(std::string  name);

        void StartTimer();

        void StopTimer();

        long GetElapsedTimeMilliseconds() const;

        double GetElapsedTimeSeconds() const;

        void PrintElapsedTime() const;
    protected:
    private:
        std::chrono::steady_clock::time_point begin_;
        std::chrono::steady_clock::time_point end_;
        long elapsed_time_;
        const std::string name_;
    };
}


#endif //BILEVEL_GAIT_GEN_TIMER_H
