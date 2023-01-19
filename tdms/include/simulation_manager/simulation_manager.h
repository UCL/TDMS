/**
 * @file simulation_manager.h
 * @brief Class that runs the physics of the TDMS executable, and produces an output.
 */
#pragma once

#include <complex>
#include <vector>

#include "arrays.h"
#include "globals.h"
#include "objects_from_infile.h"
#include "cell_coordinate.h"
#include "loop_timers.h"
#include "pstd_variables.h"
#include "fdtd_bootstrapper.h"
#include "output_matrices.h"
#include "matrix.h"

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
  ObjectsFromInfile inputs;
  LoopTimers timers;
  PSTDVariables PSTD;
  FDTDBootstrapper FDTD;
  OutputMatrices outputs;

  PreferredInterpolationMethods pim;
  SolverMethod solver_method;

  EHVec eh_vec;

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

  void extract_phasor_norms(int frequency_index, int tind, int Nt);

  void prepare_output(const mxArray *fieldsample, const mxArray* campssample);


public:
  SimulationManager(InputMatrices in_matrices, SolverMethod _solver_method, PreferredInterpolationMethods _pim);

  IJKDims n_Yee_cells() { return inputs.IJK_tot; }

  void execute();

  OutputMatrices &get_outputs() { return outputs; }
};
