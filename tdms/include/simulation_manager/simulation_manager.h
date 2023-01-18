/**
 * @file simulation_manager.h
 * @brief
 */
#pragma once

#include <complex>
#include <vector>

#include "cell_coordinate.h"
#include "fdtd_bootstrapper.h"
#include "globals.h"
#include "loop_timers.h"
#include "loop_variables.h"
#include "matrix.h"
#include "objects_from_infile.h"
#include "output_matrices.h"
#include "pstd_variables.h"

struct norm_vector {
  std::vector<std::complex<double>> E;
  std::vector<std::complex<double>> H;
};

class SimulationManager {
private:
  //parameter to control the width of the ramp when introducing the waveform in steady state mode
  double ramp_width = 4.;

  ObjectsFromInfile inputs;
  LoopTimers timers;
  PSTDVariables PSTD;
  FDTDBootstrapper FDTD;
  OutputMatrices outputs;

  PreferredInterpolationMethods pim;
  SolverMethod solver_method;

  std::vector<std::complex<double>>
          E_norm,//< Holds the E-field phasors norm at each extraction frequency
          H_norm;//< Holds the H-field phasors norm at each extraction frequency

  /**
   * @brief Prepare the simulation output
   *
   * @param fieldsample The fieldsample data from the input
   * @param campssample User-specified vertices at which to extract phasors, from the input
   */
  void prepare_output(const mxArray *fieldsample, const mxArray* campssample);


  /*Implements a linear ramp which has the properties:

  ramp(t) = 1 if t > rampwidth*period
  = t/(rampwidth*period) otherwise

  t - the current time at which to evaluate the ramp
  period - the period of the monochormatic sinusoidal excitation
  rampwidth - the fraction or number of periods with of the ramp

*/

  /**
   * @brief Impliments the linear ramp function ramp(t)
   *
   * ramp(t) = t/(ramp_width * period)   if t <= rampwidth * period,
   *         = 1                         otherwise.
   *
   * @param t Argument to ramp(t)
   * @return double ramp(t)
   */
  double linear_ramp(double t);

  /**
   * @brief Extract the E,H phasor norms at the given frequency index
   *
   * @param frequency_index Index of the frequency to extract at
   * @param tind The current timestep
   * @param Nt Number of timesteps in a sinusoidal period
   */
  void extract_phasor_norms(int frequency_index, int tind, int Nt);

public:
  SimulationManager(InputMatrices in_matrices, SolverMethod _solver_method,
                    PreferredInterpolationMethods _pim);

  IJKDims n_Yee_cells() { return inputs.IJK_tot; }

  void execute_simulation();

  OutputMatrices post_loop_processing();
};


// SimulationManager(in_matrices, solver_method, pim)
// SimulationManager.execute_simulation()
// OutputMatrices = SimulationManager.post_loop_processing()
