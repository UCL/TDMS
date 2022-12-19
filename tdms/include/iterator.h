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

/**
 * @brief Executes the main simulation.
 * Used to be the MATLAB mexFunction
 */
void execute_simulation(int, mxArray **, int, const mxArray **, SolverMethod solver_method,
                        PreferredInterpolationMethods preferred_interpolation_methods);

void extractPhasorsPlane( double **iwave_lEx_Rbs, double **iwave_lEx_Ibs, double **iwave_lEy_Rbs, double **iwave_lEy_Ibs,
			  double **iwave_lHx_Rbs, double **iwave_lHx_Ibs, double **iwave_lHy_Rbs, double **iwave_lHy_Ibs,
        ElectricSplitField &E, MagneticSplitField &H,
			  int I_tot, int J_tot, int K1, int n, double omega, double dt, int Nt);

void initialiseDouble3DArray(double ***inArray, int i_lim, int j_lim, int k_lim);

void initialiseDouble2DArray(double **inArray, int i_lim, int j_lim);

double linearRamp(double t, double period, double rampwidth);

bool is_dispersive(unsigned char ***materials,double *gamma, double dt, int I_tot, int J_tot, int K_tot);

void extractPhasorENorm(std::complex<double> *Enorm, double ft, int n, double omega, double dt, int Nt);

void extractPhasorHNorm(std::complex<double> *Hnorm, double ft, int n, double omega, double dt, int Nt);

void normaliseVertices( double **EHr, double **EHi, ComplexAmplitudeSample &campssample, std::complex<double> Enorm , std::complex<double> Hnorm );

void update_EH(double **EHr, double **EHi, int vindex, int idx, std::complex<double> &phase_term, double &value);

void extractPhasorsVertices(double **EHr, double **EHi, ElectricSplitField &E, MagneticSplitField &H,
                            ComplexAmplitudeSample &campssample, int n, double omega, double dt, int Nt,
                            int dimension,int J_tot,int intmethod );
