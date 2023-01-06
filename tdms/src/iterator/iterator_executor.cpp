#include "iterator_executor.h"

#include <ctime>
#include <spdlog/spdlog.h>

#include "globals.h"

using tdms_math_constants::DCPI;

void Iterator_Executor::log_current_iteration_and_max_fields() {
  // if it has been long enough since we last printed the current iteration number
  if ((((double) time(NULL)) - last_logged_iteration_time) > seconds_between_logs) {

    double maxfield = std::max(E_s.largest_field_value(), H_s.largest_field_value());

    spdlog::info("Iterating: tind = {0:d}, maxfield = {1:e}", tind, maxfield);
    last_logged_iteration_time = double(time(NULL));
  }
}

void Iterator_Executor::prepare_phasor_convergence_proceedure() {
  // First we set dt so that an integer number of time periods fits within a sinusoidal period
  double Nsteps_tmp = 0.0;
  double dt_old;
  if (params.source_mode == SourceMode::steadystate) {
    dt_old = params.dt;
    Nsteps_tmp = ceil(2. * DCPI / params.omega_an / params.dt * 3);
    params.dt = 2. * DCPI / params.omega_an * 3 / Nsteps_tmp;

    // in a complete run, display to the user that dt had to be changed
    if (params.run_mode == RunMode::complete) {
      spdlog::info("Changing dt from {0:.10e} to {1:.10e}", dt_old, params.dt);
    }
  }

  // set Nsteps based on previous calculation, if appropriate
  Nsteps = (int) lround(Nsteps_tmp);

  // params.Nt should be an integer number of Nsteps when operating in steady-state
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete) {
    if (params.Nt / Nsteps * Nsteps != params.Nt) {
      unsigned int tmp_Nt = params.Nt;//< Temporary storage for the log
      params.Nt = params.Nt / Nsteps * Nsteps;
      spdlog::info("Changing the value of Nt from {0:d} to {1:d} for correct phasor extraction",
                   tmp_Nt, params.Nt);
    }
    // display Nsteps
    spdlog::info("Using Nsteps = {0:d}", Nsteps);
  }

}
