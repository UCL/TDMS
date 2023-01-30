#include "simulation_manager/pstd_variables.h"

#include <omp.h>

#include "numerical_derivative.h"

void PSTDVariables::set_using_dimensions(const IJKDimensions &IJK_tot) {
  int max_IJK = IJK_tot.max(), n_threads = omp_get_max_threads();

  ca.allocate(n_threads, max_IJK + 1);
  cb.allocate(n_threads, max_IJK + 1);

  // deduce the number of coefficients in the derivative-shift operator
  N_ex = IJK_tot.i;
  N_ey = IJK_tot.j;
  N_ez = IJK_tot.k;
  N_hx = IJK_tot.i + 1;
  N_hy = IJK_tot.j + 1;
  N_hz = IJK_tot.k + 1;

  // flag allocation of memory
  variables_allocated = true;
  // allocate suitable memory for the derivative-shift operator coefficients
  dk_ex = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_ex));
  dk_ey = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_ey));
  dk_ez = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_ez));
  dk_hx = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_hx));
  dk_hy = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_hy));
  dk_hz = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_hz));

  // initialise the derivative-shift operators for each field component
  init_diff_shift_op(-0.5, dk_ex, N_ex);
  init_diff_shift_op(-0.5, dk_ey, N_ey);
  init_diff_shift_op(-0.5, dk_ez, N_ez);
  init_diff_shift_op(0.5, dk_hx, N_hx);
  init_diff_shift_op(0.5, dk_hy, N_hy);
  init_diff_shift_op(0.5, dk_hz, N_hz);
}

PSTDVariables::~PSTDVariables() {
  if (variables_allocated) {
    fftw_free(dk_ey);
    fftw_free(dk_ez);
    fftw_free(dk_hx);
    fftw_free(dk_hy);
    fftw_free(dk_hz);
  }
}
