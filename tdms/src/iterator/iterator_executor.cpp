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

void Iterator_Executor::optimise_loops_if_possible(double non_zero_tol) {
  bool ksource_nz[4];//< "k source non-zero"
  for (int icomp = 0; icomp < 4; icomp++) { ksource_nz[icomp] = false; }

  if (J_tot == 0) {
    for (int icomp = 0; icomp < 4; icomp++)
      for (int ki = 0; ki < (I_tot + 1); ki++) {
        ksource_nz[icomp] = ksource_nz[icomp] ||
                            (fabs(Ksource.imag[0][ki - (I0.index)][icomp]) > non_zero_tol) ||
                            (fabs(Ksource.real[0][ki - (I0.index)][icomp]) > non_zero_tol);
      }
  }
  /* We now know the following information:
    Ey and Hx receive an input from the source condition only if ksource_nz[2] or ksource_nz[1]
    are non-zero.
    Ex and Hy receive an input from the source condition only if ksource_nz[3] or ksource_nz[0]
    are non-zero.

    Thus we can set the variables J_tot_p1_bound, J_tot_bound accordingly for the 3D, 2D-TE, and 2D-TM simulations.
  */

  J_tot_bound = J_tot;
  J_tot_p1_bound = J_tot + 1;
  // If in a 2D simulation, adjust the values accordingly
  if (J_tot == 0) {
    // TE case
    if (ksource_nz[2] || ksource_nz[1] || params.eyi_present) {
      J_tot_bound = 1;
    } else {
      J_tot_bound = 0;
    }
    // TM case
    if (ksource_nz[3] || ksource_nz[0] || params.exi_present) {
      J_tot_p1_bound = 1;
    } else {
      J_tot_p1_bound = 0;
    }
  }
}
