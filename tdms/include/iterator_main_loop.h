#pragma once

#include <ctime>
#include <stdexcept>

#include <spdlog/spdlog.h>

#include "field.h"
#include "simulation_parameters.h"
#include "timer.h"

enum class IterationTimers { MAIN, INTERNAL };

class IteratorMainLoop {
private:
  Timer main_loop_timer;//< Timer for the main iteration loops
  Timer internal_timer;//< Timer for tasks that are internal to one iteration
  double last_logged_iteration_time = double(time(NULL));//< The absolute time at which we last logged the iteration number
  double seconds_between_logs = 1.;//< The time in seconds to wait between printing the current iteration number to the log
  // Error to throw when we try to use a timer that doesn't exist
  void unrecognised_timer() { throw std::runtime_error("Didn't recognise that timer!"); }

public:
  int dft_counter = 0;
  int Nsteps = 0;

  void start_timer(IterationTimers timer) {
    switch (timer) {
        case IterationTimers::MAIN:
            main_loop_timer.start();
            break;
        case IterationTimers::INTERNAL:
            internal_timer.start();
            break;
        default:
            unrecognised_timer();
            break;
    }
  }
  void end_timer(IterationTimers timer) {
    switch (timer) {
        case IterationTimers::MAIN:
            main_loop_timer.end();
            break;
        case IterationTimers::INTERNAL:
            internal_timer.end();
            break;
        default:
            unrecognised_timer();
            break;
    }
  }
  void click_timer(IterationTimers timer) {
    switch (timer) {
        case IterationTimers::MAIN:
            main_loop_timer.click();
            break;
        case IterationTimers::INTERNAL:
            internal_timer.click();
            break;
        default:
            unrecognised_timer();
            break;
    }
  }
  double time_ellapsed_by(IterationTimers timer) {
    switch (timer) {
        case IterationTimers::MAIN:
            return main_loop_timer.delta_seconds();
            break;
        case IterationTimers::INTERNAL:
            return internal_timer.delta_seconds();
            break;
        default:
            unrecognised_timer();
            return -1.;
            break;
    }
  }

  /**
   * @brief Print the current iteration number and max-field information, if it has been long enough since we last printed this information.
   *
   * @param E_split,H_split Split-fields currently being iterated
   * @param tind Current iteration number
   */
  void log_update(ElectricSplitField &E_split, MagneticSplitField &H_split, int tind);

};
