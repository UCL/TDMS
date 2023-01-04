#pragma once

#include <stdexcept>

#include "simulation_parameters.h"
#include "timer.h"

enum class IterationTimers { MAIN, INTERNAL };

class IteratorMainLoop {
private:
  Timer main_loop_timer, internal_timer;

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
};
