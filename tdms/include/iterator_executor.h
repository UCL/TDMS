/**
 * @file iterator_executor.h
 * @brief Class that runs the TDMS simulation and returns the desired outputs
 */
#pragma once;

#include "iterator_timers.h"
#include "iterator_loop_variables.h"

/**
 * @brief This class will run the main TDMS simulation, and then handle how the output is passed back to main().
 *
 * Variables needed for the main loop are handled by inheriting from Iterator_LoopVariables.
 *
 * Timing is controlled by inheritence from Iterator_Logger
 *
 */
class Iterator_Executor : public Iterator_LoopVariables, Iterator_Timers {
private:
  double last_logged_iteration_time =
          double(time(NULL));//< The absolute time at which we last logged the iteration number
  double seconds_between_logs =
          1.;//< The time in seconds to wait between printing the current iteration number to the log

public:
  int dft_counter =
          0;     //< Number of DFTs that have been performed since last phasor convergence check
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

  // Use the Iterator_LoopVariables constructor - Iterator_Timers default initaliser should be sufficient for remaining members
  Iterator_Executor(InputMatrices matrices_from_input_file, SolverMethod _solver_method, PreferredInterpolationMethods interpolation_method) : Iterator_LoopVariables(matrices_from_input_file, _solver_method, interpolation_method), Iterator_Timers() { spdlog::info("Building Iterator_Executor object..."); }

  /**
   * @brief Sets time_E and time_H according to the current iteration and timestep dt
   *
   * @param dt The simulation timestep
   */
  void update_field_times_to_current_iteration(double dt) {
    time_E = ((double) (tind + 1)) * dt;
    time_H = time_E - dt / 2.;
  }

  /**
   * @brief Print the current iteration number and max-field information, if it has been long enough since we last printed this information.
   */
  void log_current_iteration_and_max_fields();

  /**
     * @brief Sets up the parameters that will be required when checking phasor convergence.
     *
     * Nsteps - the number of DFTs we are allowed to perform before we check for convergence again
     *
     *
     */
  void prepare_phasor_convergence_proceedure();
  //void run_main_loop();
};
