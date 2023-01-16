#pragma once

#include <fftw3.h>

#include "arrays.h"
#include "cell_coordinate.h"

/* PSTD EXCLUSIVE VARIABLES (REQUIRES MEMORY MANAGEMENT)
  These variables are not used when using FDTD method */
class PSTDVariables {
private:
  // flags whether the dk-variables have been malloc'd or not
  bool variables_allocated = false;

public:
  CCoefficientMatrix ca_vec, cb_vec, cc_vec;

  // The number of complex fourier coefficients in the derivative-shift operator, for each field component ( N_e_x = # coeffs for Ex, for example)
  int N_e_x, N_e_y, N_e_z, N_h_x, N_h_y, N_h_z;
  // The coefficients of the PSTD derivative shift operator, for each field component ( dk_e_x = coeffs for Ex, for example )
  fftw_complex *dk_e_x, *dk_e_y, *dk_e_z, *dk_h_x, *dk_h_y, *dk_h_z;

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
