/**
 * @file loop_timers.h
 * @brief Declares the class that handles the timing of the main loop and
 * subprocesses of the main loop.
 */
#pragma once

#include <ctime>
#include <stdexcept>

#include "timer.h"

// There are two timers concerning the main loop: the internal timer that is
// tracking individual tasks, and the timer that is tracking how long the main
// iteration has actually taken.
enum class TimersTrackingLoop { MAIN, INTERNAL };

/**
 * @brief Class that handles the timers that monitor the performance of the
 * program during the main loop, and providing information to be printed to the
 * logger/debugger.
 */
class LoopTimers {
private:
  Timer main_loop_timer;//< Timer for the main iteration loops
  Timer internal_timer; //< Timer for tasks that are internal to one iteration
  // Error to throw when we try to use a timer that doesn't exist
  void unrecognised_timer() {
    throw std::runtime_error("Didn't recognise that timer!");
  }

public:
  void start_timer(TimersTrackingLoop timer) {
    switch (timer) {
      case TimersTrackingLoop::MAIN:
        main_loop_timer.start();
        break;
      case TimersTrackingLoop::INTERNAL:
        internal_timer.start();
        break;
      default:
        unrecognised_timer();
        break;
    }
  }
  void end_timer(TimersTrackingLoop timer) {
    switch (timer) {
      case TimersTrackingLoop::MAIN:
        main_loop_timer.end();
        break;
      case TimersTrackingLoop::INTERNAL:
        internal_timer.end();
        break;
      default:
        unrecognised_timer();
        break;
    }
  }
  void click_timer(TimersTrackingLoop timer) {
    switch (timer) {
      case TimersTrackingLoop::MAIN:
        main_loop_timer.click();
        break;
      case TimersTrackingLoop::INTERNAL:
        internal_timer.click();
        break;
      default:
        unrecognised_timer();
        break;
    }
  }
  double time_ellapsed_by(TimersTrackingLoop timer) {
    switch (timer) {
      case TimersTrackingLoop::MAIN:
        return main_loop_timer.delta_seconds();
        break;
      case TimersTrackingLoop::INTERNAL:
        return internal_timer.delta_seconds();
        break;
      default:
        unrecognised_timer();
        // return garbage just to make the compiler happy - the statement above
        // will error so we'll never get here anyway
        return -1.;
        break;
    }
  }
};
