/**
 * @file timer.h
 * @brief Class wrapper for timing with OpenMP's wall timing.
 */
#pragma once

/**
 * @brief Stopwatch class.
 */
class Timer{

    double start_time; //< start time in seconds
    double end_time;   //< end time in seconds

public:
    /** Starts the stopwatch */
    void start();
    /** Stops the stopwatch */
    void end();
    /** Log the difference in time and reset the timer */
    void click();
    /** Time difference */
    double delta_seconds() const;
};
