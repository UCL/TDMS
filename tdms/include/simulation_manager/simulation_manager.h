/**
 * @file simulation_manager.h
 * @brief
 */
#pragma once

#include "arrays.h"
#include "globals.h"
#include "objects_from_infile.h"
#include "cell_coordinate.h"
#include "loop_timers.h"
#include "pstd_variables.h"
#include "fdtd_bootstrapper.h"
#include "output_matrices.h"
#include "matrix.h"

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

  void prepare_output(const mxArray *fieldsample, const mxArray* campssample);

public:
  SimulationManager(InputMatrices in_matrices, SolverMethod _solver_method, PreferredInterpolationMethods _pim);

  IJKDims n_Yee_cells() { return inputs.IJK_tot; }

};
