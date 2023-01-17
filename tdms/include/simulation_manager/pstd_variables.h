/**
 * @file pstd_variables.h
 * @brief Declares the PSTDVariables class, which handles variables exclusively used in the PSTD solver method
 */
#pragma once

#include <fftw3.h>

#include "arrays.h"
#include "cell_coordinate.h"

/**
 * @brief Handles allocation and tear-down of memory/variables required when running a PSTD simulation. If running an FDTD simulation, none of the member variables of this class are needed.
 */
class PSTDVariables {
private:
  // flags whether the dk-variables have been malloc'd or not
  bool variables_allocated = false;

public:
  CCoefficientMatrix ca, cb;//< TODO

  // The number of complex fourier coefficients in the derivative-shift operator, for each field component ( N_ex = # coeffs for Ex, for example)
  int N_ex, N_ey, N_ez, N_hx, N_hy, N_hz;
  // The coefficients of the PSTD derivative shift operator, for each field component ( dk_ex = coeffs for Ex, for example )
  fftw_complex *dk_ex, *dk_ey, *dk_ez, *dk_hx, *dk_hy, *dk_hz;

  PSTDVariables() = default;
  PSTDVariables(IJKDims IJK_tot) { set_using_dimensions(IJK_tot); }

  /**
       * @brief Allocate memory for PSTD method, for a simulation with the provided number of Yee cells in each dimension.
       *
       * @param IJK_tot Triple containing the number of Yee cells in the I,J,K directions
       */
  void set_using_dimensions(IJKDims IJK_tot);

  ~PSTDVariables();
};
