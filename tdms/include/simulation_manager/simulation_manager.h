/**
 * @file simulation_manager.h
 * @brief Class that runs the physics of the TDMS executable, and produces an
 * output. Implementation is split between simulation_manager.cpp and
 * execute_simulation.cpp.
 */
#pragma once

#include <complex>
#include <string>
#include <vector>

#include "arrays.h"
#include "cell_coordinate.h"
#include "globals.h"
#include "input_flags.h"
#include "input_matrices.h"
#include "output_matrices/output_matrices.h"
#include "simulation_manager/fdtd_bootstrapper.h"
#include "simulation_manager/loop_timers.h"
#include "simulation_manager/loop_variables.h"
#include "simulation_manager/objects_from_infile.h"
#include "simulation_manager/pstd_variables.h"

// Whether or not to time execution of loop subtasks
#define TIME_EXEC false
// Whether or not to time the main loop execution
#define TIME_MAIN_LOOP true
// Threshold used to terminate the steady state iterations
#define TOL 1e-6

/**
 * @brief Manages the physics of TDMS and the simulation loop itself.
 *
 * Takes in the contents of the input file and produces the output in the
 * outputs member, so that it can be written to the specified output file.
 *
 * Attributes are private to avoid external alterations to their content and
 * potential reassignment of MATLAB pointers, resulting in leaking memory. The
 * outputs can be accessed through a fetch-method within main.cpp if necessary.
 */
class SimulationManager {
private:
  /*! The solver method to use in this simulation */
  tdms_flags::SolverMethod solver_method;
  /*! The interpolation methods to use in this simulation */
  tdms_flags::InterpolationMethod i_method;

  /*! The input objects that are generated from an input file */
  ObjectsFromInfile inputs;

  LoopTimers timers; //!< Timers for tracking the execution of the simulation
  PSTDVariables PSTD;//!< PSTD-solver-specific variables

  /*! Output object that will contain the results of this simulation, given the
   * input file */
  OutputMatrices outputs;

  EHVec eh_vec;//!< TODO

  /*! Width of the ramp when introducing the waveform in steady state mode */
  double ramp_width = 4.;
  /**
   * @brief Evaluates the linear ramp function r(t)
   *
   * r(t) = 1                           if t > ramp_width * period
   *      = t / (ramp_width * period)   if t <= ramp_width * period
   *
   * Period = 1. / (omega_an / (2 * pi)) = 2 * pi / omega_an
   *
   * NB Note that r(t) = min(1, t / (ramp_width * period))
   *
   * @param t Time to evaluate function at
   * @return double
   */
  double linear_ramp(double t) {
    double period = 2 * tdms_math_constants::DCPI / inputs.params.omega_an;
    return std::min(1., t / (ramp_width * period));
  }

  /*! Holds the {E,H}-field phasors norm at each extraction frequency */
  std::vector<std::complex<double>> E_norm;
  /*! @copydoc E_norm */
  std::vector<std::complex<double>> H_norm;
  /**
   * @brief Extracts the phasors norms at the given frequency (index)
   *
   * @param frequency_index The index of the frequency (in inputs.f_ex_vec) to
   * extract at
   * @param tind The current timestep
   * @param Nt The number of timesteps in a sinusoidal period
   */
  void extract_phasor_norms(int frequency_index, unsigned int tind, int Nt);

  /**
   * @brief Creates MATLAB memory blocks that will be iteratively updated in the
   * execute() method, and assigns array sizes for writing purposes later.
   *
   * @param fieldsample The fieldsample input from the input file
   * @param campssample The complex amplitude sample input from the input file
   */
  void prepare_output(const mxArray *fieldsample, const mxArray *campssample);

  /**
   * @brief Update electric-split field components and current densities AT A
   * PARTICULAR CELL in steady-state after an E-field timestep has been
   * performed.
   *
   * @param time_H The time the magnetic field is currently sitting at.
   * @param parallel The axis to which the plane of the Source term is parallel.
   * EG X = Isource, Y = Jsource, etc
   * @param C_axis Whether we are updating terms along the C-axis (true) or
   * B-axis (false) of the Source plane. See Source doc for axis information.
   * @param zero_plane Whether we are updating terms on the 0-plane (true) or
   * the 1-plane (false). EG I0 would have this input as true, whereas J1 as
   * false.
   * @param is_conductive Whether the medium is conductive (so J_c needs to be
   * updated)
   * @param cell_b,cell_c The coordinates (cell_a, cell_b, cell_c) of the Yee
   * cell in which we are updating the field. See the Source doc for notation
   * information.
   * @param array_ind The index of the various material property arrays that
   * correspond to this particular Yee cell, and thus update equation.
   * @param J_c The current density in the conductive medium
   * @param J_s The current density in the dispersive medium
   */
  void E_source_update_steadystate(double time_H, AxialDirection parallel,
                                   bool C_axis, bool zero_plane,
                                   bool is_conductive, int cell_b, int cell_c,
                                   int array_ind, CurrentDensitySplitField &J_c,
                                   CurrentDensitySplitField &J_s);
  /**
   * @brief Update magnetic-split field components AT A
   * PARTICULAR CELL in steady-state after an H-field timestep has been
   * performed.
   *
   * @param time_E The time the electric field is currently sitting at
   * @param parallel The axis to which the plane of the Source term is parallel.
   * EG X = Isource, Y = Jsource, etc
   * @param C_axis Whether we are updating terms along the C-axis (true) or
   * B-axis (false) of the Source plane. See Source doc for axis information.
   * @param zero_plane Whether we are updating terms on the 0-plane (true) or
   * the 1-plane (false). EG I0 would have this input as true, whereas J1 as
   * false.
   * @param array_ind The index of the various material property arrays that
   * correspond to this particular Yee cell, and thus update equation.
   * @param cell_b,cell_c The coordinates (cell_a, cell_b, cell_c) of the Yee
   * cell in which we are updating the field. See the Source doc for notation
   * information.
   */
  void H_source_update_steadystate(double time_E, AxialDirection parallel,
                                   bool zero_plane, bool C_axis, int array_ind,
                                   int cell_b, int cell_c);
  /**
   * @brief [E-FIELD UPDATES] Performs updates to the electric-split field
   * components and current density fields after an E-field timestep has been
   * performed, in accordance with the I,J, and K-source terms.
   *
   * @param time_H The time the magnetic field is currently sitting at.
   * @param is_conductive Whether the medium is conductive (so J_c needs to be
   * updated)
   * @param J_c The current density in the conductive medium
   * @param J_s The current density in the dispersive medium
   */
  void E_source_update_all_steadystate(double time_H, bool is_conductive,
                                       CurrentDensitySplitField &J_c,
                                       CurrentDensitySplitField &J_s);
  /**
   * @brief [H-FIELD UPDATES] Performs updates to the magnetic-split field
   * components after an H-field timestep has been performed, in accordance with
   * the I,J, and K-source terms.
   *
   * @param time_E The time the electric field is currently sitting at.
   */
  void H_source_update_all_steadystate(double time_E);

  /*! @copydoc E_source_update_all_steadystate */
  void E_Isource_update_steadystate(double time_H, bool is_conductive,
                                    CurrentDensitySplitField &J_c,
                                    CurrentDensitySplitField &J_s);
  /*! @copydoc E_source_update_all_steadystate */
  void E_Jsource_update_steadystate(double time_H, bool is_conductive,
                                    CurrentDensitySplitField &J_c,
                                    CurrentDensitySplitField &J_s);
  /*! @copydoc E_source_update_all_steadystate */
  void E_Ksource_update_steadystate(double time_H, bool is_conductive,
                                    CurrentDensitySplitField &J_c,
                                    CurrentDensitySplitField &J_s);
  /*! @copydoc H_source_update_all_steadystate */
  void H_Isource_update_steadystate(double time_E);
  /*! @copydoc H_source_update_all_steadystate */
  void H_Jsource_update_steadystate(double time_E);
  /*! @copydoc H_source_update_all_steadystate */
  void H_Ksource_update_steadystate(double time_E);

  /**
   * @brief [E-FIELD UPDATES] Performs the updates to the electric-split field
   * components nad current density fields after an E-field timestep has been
   * performed, in accordance with the source terms.
   *
   * @param time_H The time the magnetic field is currently sitting at.
   * @param is_conductive Whether the medium is conductive (so J_c needs to be
   * updated)
   * @param J_c The current density in the conductive medium
   * @param J_s The current density in the dispersive medium
   */
  void update_source_terms_pulsed(double time_H, bool is_conductive,
                                  CurrentDensitySplitField &J_c,
                                  CurrentDensitySplitField &J_s);
  /**
   * @brief [H-FIELD UPDATES] Performs the updates to the magnetic-split field
   * components and current density fields after an H-field timestep has been
   * performed, in accordance with the source terms.
   *
   * @param time_E The time the electric field is currently sitting at.
   * @param tind The current iteration number
   */
  void update_source_terms_pulsed(double time_E, int tind);

  /* execute() subfunctions to break up main loop */

  /**
   * @brief Checks whether the phasors have converged in a steady-state
   * simulation.
   *
   * E and H fields are zero'd in this method if convergence is reached. E_copy
   * is updated to hold the converged values in this case.
   *
   * @param[inout] dft_counter The number of DFTs that have been performed since
   * we began checking for convergence
   * @param[inout] E_copy The field array that stores the phasors from the
   * previous iteration.
   * @return true If the phasors have converged
   * @return false Phasors have not converged
   */
  bool check_phasor_convergence(int &dft_counter, ElectricField &E_copy);
  /**
   * @brief Extracts the phasors in the volume, on the user-defined surface, and
   * at user-defined vertices.
   *
   * If the user has not requested one or more of these extraction locations,
   * the corresponding steps are skipped. If the RunMode is not complete, then
   * this step is always bypassed.
   *
   * @param dft_counter The number of DFTs that have been performed since we
   * began checking for convergence
   * @param tind The current iteration number
   */
  void extract_phasors(int &dft_counter, unsigned int tind);
  /**
   * @brief Computes the detector function.
   *
   * Presumably this is the functional form of the fields/phasors at the
   * detector positions.
   *
   * @param tind The current iteraton number
   * @param lv Variables required from the main loop
   */
  void compute_detector_functions(unsigned int tind, LoopVariables &lv);
  /**
   * @brief Begin a new acquisition period: zero the angular-norms and field
   * normalisation factors.
   *
   * Each time a new acquisition period of harmonic illumination begins, all
   * complex amplitudes (volume, surface etc.) are set back to 0. This is
   * because the discrete Fourier transforms used to acquire these complex
   * amplitudes starts again.
   *
   * In particular, the returned complex amplitudes will have been acquired
   * during a single acquisition period of harmonic illumination. Note that the
   * acquisition period is actually three periods of the harmonic waves'
   * fundamental period.
   *
   * @param tind The current iteration number
   */
  void new_acquisition_period(unsigned int tind);
  /**
   * @brief Run tasks at the end of an iteration, in preparation for the next.
   *
   * These tasks are:
   * - Write to the log with the field residual if it has been long enough since
   * the last write
   * - Report on possible convergence failure (and set outputs accordingly)
   * - Export fields if at a suitable iteration number
   *
   * @param[inout] time_of_last_log Time since we last wrote to the log file. If
   * we write to the log file whilst running this method, this value is updated
   * to the time of writing.
   * @param[in] tind The current iteration number.
   * @param[inout] E_copy The object that is storing the phasors from the
   * previous iteration, for use in convergence checking.
   */
  void end_of_iteration_steps(double &time_of_last_log, unsigned int tind,
                              ElectricField &E_copy);

  /* execute() subfunctions that time-propagate fields */

  void update_Exy(LoopVariables &lv);
  void update_Exz(LoopVariables &lv);

public:
  SimulationManager(InputMatrices in_matrices, const InputFlags &in_flags);

  /** @brief Fetch the number of Yee cells in each dimension */
  IJKDimensions n_Yee_cells() { return inputs.IJK_tot; }

  /** @brief Run the time-stepping algorithm given the current inputs. */
  void execute();

  /**
   * @brief Perform additional (mostly conditional) computations on the outputs,
   * depending on the simulation inputs and run specifications.
   *
   * This should only be run AFTER a successful run of the execute() method.
   */
  void post_loop_processing();

  /**
   * @brief Write the outputs to the file provided. Wrapper for
   * outputs.save_outputs
   *
   * @param output_file The filename to write the outputs to
   * @param compressed_format If true, write compressed output (do not write
   * facets and vertices)
   */
  void write_outputs_to_file(std::string output_file,
                             bool compressed_format = false) {
    outputs.save_outputs(output_file, compressed_format);
  }
};
