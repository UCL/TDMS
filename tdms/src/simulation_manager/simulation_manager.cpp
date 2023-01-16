#include "simulation_manager.h"

SimulationManager::SimulationManager(InputMatrices inputs, SolverMethod _solver_method)
    : infile_data(inputs, _solver_method) {}
