/**
 * @file simulation_manager.h
 * @brief
 */
#pragma once

#include "globals.h"
#include "objects_from_infile.h"
#include "cell_coordinate.h"
#include "loop_timers.h"
#include "pstd_variables.h"
#include "fdtd_bootstrapper.h"

class SimulationManager {
private:
  ObjectsFromInfile infile_data;
  LoopTimers timers;
  PSTDVariables pstd;
  FDTDBootstrapper fdtd_boots;

public:
  SimulationManager(InputMatrices inputs, SolverMethod _solver_method);

  IJKTotal n_Yee_cells() { return infile_data.IJK_tot; }
};
