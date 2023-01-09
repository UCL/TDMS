/**
 * @file iterator_executor.h
 * @brief Class that runs the TDMS simulation and returns the desired outputs
 */
#pragma once

#include <spdlog/spdlog.h>

#include "iterator_loop_variables.h"
#include "iterator_timers.h"

//whether of not to time execution
#define TIME_EXEC false
// time main loop
#define TIME_MAIN_LOOP true
//threshold used to terminate the steady state iterations
#define TOL 1e-6
//parameter to control the PreferredInterpolationMethodswith of the ramp when introducing the waveform in steady state mode
#define ramp_width 4.

/**
 * @brief This class will run the main TDMS simulation, and then handle how the output is passed back to main().
 *
 * Variables needed for the main loop are handled by inheriting from Iterator_LoopVariables.
 *
 * Timing is controlled by inheritence from Iterator_Logger
 *
 */
class Iterator_Executor : public Iterator_LoopVariables, public Iterator_Timers {
private:
  double last_logged_iteration_time =
          double(time(NULL));//< The absolute time at which we last logged the iteration number
  double seconds_between_logs =
          1.;//< The time in seconds to wait between printing the current iteration number to the log and the previous print

  /* DECLARATIONS FOR SUBFUNCTIONS OF THE MAIN LOOP */

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
   * @brief Implements a linear ramp function.
   *
   * This is the function
   * ramp(t) =  1                           if t > ramp_width * period,
   *         =  t / (ramp_width * period)   otherwise.
   *
   * ramp_width is defined in the header file iterator_executor.h.
   *
   * @param t Time at which to evaluate the ramp
   * @param period Period of the monochromatic sinusoidal excitation
   * @return double ramp(t), defined above
   */
  double linear_ramp(double t, double period) {
    if (t > period * ramp_width) {
      return 1.;
    } else {
      return t / (period * ramp_width);
    }
  }

  /**
   * @brief Determine whether the main loop should be terminated due to phasor convergence being achieved.
   *
   * If convergence has not been achieved, setup the next iteration accordingly before exiting this method.
   *
   * @return true Phasor convergence achieved - terminate main loop immediately.
   * @return false Phasor convergence has not been achieved - setup the next iteration before returning.
   */
  bool phasors_have_converged();

  /**
   * @brief Extract the phasor norm of the {E,H}-field at the given frequency to {E,H}_norms
   *
   * @param frequency_index Index of the frequency to extract at
   * @param Nt Number of timesteps in a single sinusoidal period ( ? -> CONFIRM W/ PETER)
   */
  void extract_E_phasor_norm(int frequency_index, int Nt);
  /*! @copydoc extract_E_phasor_norm */
  void extract_H_phasor_norm(int frequency_index, int Nt);
  /**
   * @brief Extract the phasor norm, across all frequencies, of both the E and H field and update the corresponding {E,H}_norm arrays.
   *
   * @param Nt Number of timesteps in a single sinusoidal period ( ? -> CONFIRM W/ PETER)
   */
  void extract_phasor_norms_at_all_frequencies(int Nt) {
    for (int frequency_index = 0; frequency_index < f_ex_vec.size(); frequency_index++) {
      extract_E_phasor_norm(frequency_index, Nt);
      extract_H_phasor_norm(frequency_index, Nt);
    }
  }

    // TODO: Docstring!
  void extract_phasors_in_plane();

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

    We note that the extract_E_phasor_norm uses (tind+1)*dt and extract_H_phasor_norm uses
    (tind+1/2)*dt to perform the update equation in the DFT. This seems incorrect
    at first but we note that they take the terms fte and fth as inputs respectively.
    When one notes that fte is calculated using time_E and fth using time_H we see
    that this indexing is correct, ie, time_E = (tind+1)*dt and time_H = (tind+1/2)*dt.
  */
  double time_E = 0., time_H = 0.;

  // Use the Iterator_LoopVariables constructor - Iterator_Timers default initaliser should be sufficient for remaining members
  Iterator_Executor(InputMatrices matrices_from_input_file, SolverMethod _solver_method,
                    PreferredInterpolationMethods interpolation_method)
      : Iterator_LoopVariables(matrices_from_input_file, _solver_method, interpolation_method),
        Iterator_Timers() {}

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

  /**
   * @brief Performs an optimisation step in a 2D simulation (J_tot==0) when we have either TE or TM, but not both.
   *
   * In the J_tot==0 2D version; the 'TE' case involves components Ey, Hx and Hz, whilst the 'TM' case involves components Ex, Ez and Hy.
   *
   * The idea is to use an alternative upper limit to the loop over j when we have J_tot==0. As it stands, we are planning the loops listed below. We use the syntax (k k_min k_max : j j_min j_max : i i_min i_max) to represent the nested for loop:
   *    for( k = k_min; k < k_max; k++) {
   *        for( j = j_min; j < j_max; j++) {
   *            for( i = i_min; i < i_max; i++) {
   *            }
   *        }
   *    }
   *
   * Exy: Not involved in 2D, (k 0 K_tot+1 : j 1 J_tot : i 0 I_tot)
   * Exz: TM, (k 1 K_tot : j 0 J_tot+1 : i 0 I_tot)
   * Eyx: TE, (k 0 K_tot+1 : j 0 max(J_tot, 1) : i 1 I_tot)
   * Eyz: TE, (k 1 K_tot : j 0 max(J_tot, 1) : i 0 I_tot+1)
   * Ezx: TM, (k 0 K_tot : j 0 J_tot+1 : i 1 I_tot)
   * Ezy: Not involved in 2D, (k 0 K_tot : j 1 J_tot : i 0 I_tot+1)
   * Hxy: Not involved in 2D, (k 0 K_tot : j 0 J_tot : i 0 I_tot+1)
   * Hxz: TE, (k 0 K_tot : j 0 max(J_tot,1) : i 0 I_tot+1)
   * Hyx: TM, (k 0 K_tot : j 0 J_tot+1 : i 0 I_tot)
   * Hyz: TM, (k 0 K_tot : j 0 J_tot+1 : i 0 I_tot)
   * Hzx: TE, (k 0 K_tot+1 : j 0 max(J_tot,1) : i 0 I_tot)
   * Hzy: Not involved in 2D, (k 0 K_tot+1) : j 0 J_tot : i 0 I_tot)
   *
   * We see that in all cases, the TE update has the following loop on j:
   *    for(j=0;j<max(J_tot,1);j++)
   * whilst the TM case has:
   *    for(j=0;j<(J_tot+1);j++)
   *
   * So we can create variables
   * int J_tot_p1_bound, J_tot_bound
   * which would take the following values
   *    3D:
   *        J_tot_p1_bound = J_tot + 1;
   *        J_tot_bound = J_tot;
   *    2D:
   *        TE:
   *            J_tot_bound = 1;
   *        Not TE:
   *            J_tot_bound = 0;
   *        TM:
   *            J_tot_p1_bound = 1;
   *        Not TM:
   *            J_tot_p1_bound = 0;
   *
   * @param non_zero_tol Tolerance at which doubles are considered "non-zero"
   */
  void optimise_loops_if_possible(double non_zero_tol = 1.e-15);

  /**
   * @brief Run the simulation main loop.
   *
   * That is, run the simulation part of the TDMS simulation.
   */
  void run_main_loop();
};
