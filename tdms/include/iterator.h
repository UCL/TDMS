/**
 * @file iterator.h
 * @brief The main time propogation algorithm.
 *
 * Contains the main FDTD loop as well as other functions such as phasor
 * extraction etc. Works in both pulsed and steady state mode.
 */
#pragma once

#include <complex>

#include "arrays.h"
#include "field.h"
#include "grid_labels.h"
#include "input_matrices.h"
#include "output_matrices.h"

/**
 * @brief Executes the main simulation.
 * Used to be the MATLAB mexFunction
 */
OutputMatrices execute_simulation(InputMatrices in_matrices, SolverMethod solver_method,
                                  PreferredInterpolationMethods preferred_interpolation_methods);

void extractPhasorsPlane(Matrix<std::complex<double>> &iwave_lEx_bs,
                         Matrix<std::complex<double>> &iwave_lEy_bs,
                         Matrix<std::complex<double>> &iwave_lHx_bs,
                         Matrix<std::complex<double>> &iwave_lHy_bs, ElectricSplitField &E,
                         MagneticSplitField &H, int I_tot, int J_tot, int K1, int n, double omega,
                         double dt, int Nt);

void initialiseDouble3DArray(double ***inArray, int i_lim, int j_lim, int k_lim);

void initialiseDouble2DArray(double **inArray, int i_lim, int j_lim);

double linearRamp(double t, double period, double rampwidth);

bool is_dispersive(unsigned char ***materials,double *gamma, double dt, int I_tot, int J_tot, int K_tot);

void extractPhasorENorm(std::complex<double> *Enorm, double ft, int n, double omega, double dt, int Nt);

void extractPhasorHNorm(std::complex<double> *Hnorm, double ft, int n, double omega, double dt, int Nt);
