/**
 * @file simulation_manager.h
 * @brief Class that runs the physics of the TDMS executable, and produces an output.
 */
#pragma once

#include <complex>
#include <vector>
#include <string>

#include "arrays.h"
#include "globals.h"
#include "objects_from_infile.h"
#include "cell_coordinate.h"
#include "loop_timers.h"
#include "pstd_variables.h"
#include "fdtd_bootstrapper.h"
#include "output_matrices.h"

/**
 * @brief Manages the physics of TDMS and the simulation loop itself.
 *
 * Takes in the contents of the input file and produces the output in the outputs member, so that it can be written to the specified output file.
 *
 * Attributes are private to avoid external alterations to their content and potential reassignment of MATLAB pointers, resulting in leaking memory. The outputs can be accessed through a fetch-method within main.cpp if necessary.
 */
class SimulationManager {
private:
  SimulationParameters params;//< The parameters for this simulation
  ObjectsFromInfile inputs;//< The input objects that are generated from an input file
  LoopTimers timers;//< Timers for tracking execution of the simulation
  PSTDVariables PSTD;//< PSTD-solver-specifi variables
  FDTDBootstrapper FDTD;//< FDTD bootsrapping variables
  OutputMatrices outputs;//< Output object that will contain the results of this simulation, given the input file

  PreferredInterpolationMethods pim;//< The interpolation methods to use in this simulation
  SolverMethod solver_method;//< The solver method to use in this simulation

  EHVec eh_vec;//< TODO

  double ramp_width = 4.;//< Width of the ramp when introducing the waveform in steady state mode
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

  std::vector<std::complex<double>> E_norm;//< Holds the E-field phasors norm at each extraction frequency
  std::vector<std::complex<double>> H_norm;//< Holds the H-field phasors norm at each extraction frequency
  /**
   * @brief Extracts the phasors norms at the given frequency (index)
   *
   * @param frequency_index The index of the frequency (in inputs.f_ex_vec) to extract at
   * @param tind The current timestep
   * @param Nt The number of timesteps in a sinusoidal period
   */
  void extract_phasor_norms(int frequency_index, int tind, int Nt);

  /**
   * @brief Creates MATLAB memory blocks that will be iteratively updated in the execute() method, and assigns array sizes for writting purposes later.
   *
   * @param fieldsample The fieldsample input from the input file
   * @param campssample The campssample input from the input file
   */
  void prepare_output(const mxArray *fieldsample, const mxArray* campssample);

public:
  SimulationManager(InputMatrices in_matrices, SolverMethod _solver_method, PreferredInterpolationMethods _pim);

  /** @brief Fetch the number of Yee cells in each dimension */
  IJKDims n_Yee_cells() { return inputs.IJK_tot; }

  /** @brief Run the time-stepping algorithm given the current inputs. */
  void execute();

  /**
   * @brief Perform additional (mostly conditional) computations on the outputs, depending on the simulation inputs and run specifications.
   *
   * This should only be run AFTER a successful run of the execute() method.
   */
  void post_loop_processing();

  /**
   * @brief Write the outputs to the file provided. Wrapper for outputs.save_outputs
   *
   * @param output_file The filename to write the outputs to
   * @param compressed_format If true, write compressed output (do not write facets and vertices)
   */
  void write_outputs_to_file(std::string output_file, bool compressed_format = false) {
    outputs.save_outputs(output_file, compressed_format);
  }
};
