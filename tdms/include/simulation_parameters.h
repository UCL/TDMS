/**
 * @file simulation_parameters.h
 * @brief Classes collecting parameters for the simulation.
 */
#pragma once

/**
 * @brief The full collection of simulation parameters.
 */
class SimulationParameters{

public:
    SimulationParameters();

    double       omega_an = 0.0;          //< Angular ω
    unsigned int Nsteps   = 0;            //< Number of simulation steps
    double       dt       = 0.0;          //< Time step
};
