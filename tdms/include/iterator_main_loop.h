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
  int dft_counter = 0;//< Number of DFTs that have been performed since last phasor convergence check
  int Nsteps = 0;//< Number of DFTs after which we should check for phasor convergence
  unsigned int tind = 0;//< Current iteration number

  /*The times of the E and H fields at the point where update equations are applied.
    time_H is actually the time of the H field when the E field consistency update is
    applied and vice versa. time_E > time_H below since after the E field consistency
    update the E field will have advanced one time step.

    The interpretation of time is slightly complicated in the following. In what follows
    I write (tind*dt,(tind+1/2)*dt) to mean the time at which we currently know the
    electric (tind*dt) and magnetic ( (tind+1/2)*dt ) fields.

    Times before                Operation         Times after
    (tind*dt,(tind+1/2)*dt)     Extract phasors   (tind*dt,(tind+1/2)*dt)
    (tind*dt,(tind+1/2)*dt)     E field update    ( (tind+1)*dt,(tind+1/2)*dt)
    ((tind+1)*dt,(tind+1/2)*dt) H field update    ( (tind+1)*dt,(tind+3/2)*dt)
    ((tind+1)*dt,(tind+3/2)*dt) Normalisation extraction

    We note that the extractPhasorENorm uses (tind+1)*dt and extractPhasorHNorm uses
    (tind+1/2)*dt to perform the update equation in the DFT. This seems incorrect
    at first but we note that they take the terms fte and fth as inputs respectively.
    When one notes that fte is calculated using time_E and fth using time_H we see
    that this indexing is correct, ie, time_E = (tind+1)*dt and time_H = (tind+1/2)*dt.
  */
  double time_E = 0., time_H = 0.;
  /**
   * @brief Sets time_E and time_H according to the current iteration and timestep dt
   *
   * @param dt The simulation timestep
   */
  void update_field_times_to_current_iteration(double dt) {
    time_E = ((double) (tind + 1)) * dt;
    time_H = time_E - dt / 2.;
  }

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
