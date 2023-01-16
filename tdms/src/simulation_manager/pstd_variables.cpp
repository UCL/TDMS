#include "pstd_variables.h"

#include <omp.h>

#include "numerical_derivative.h"

void PSTDVariables::set_using_dimensions(IJKTotal IJK_tot) {
  int max_IJK = IJK_tot.max_tot(), n_threads = omp_get_max_threads();

  ca_vec.allocate(n_threads, max_IJK + 1);
  cb_vec.allocate(n_threads, max_IJK + 1);
  cc_vec.allocate(n_threads, max_IJK + 1);

  // deduce the number of coefficients in the derivative-shift operator
  N_e_x = IJK_tot.I_tot();
  N_e_y = IJK_tot.J_tot();
  N_e_z = IJK_tot.K_tot();
  N_h_x = IJK_tot.I_tot() + 1;
  N_h_y = IJK_tot.J_tot() + 1;
  N_h_z = IJK_tot.K_tot() + 1;

  // flag allocation of memory
  variables_allocated = true;
  // allocate suitable memory for the derivative-shift operator coefficients
  dk_e_x = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_e_x));
  dk_e_y = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_e_y));
  dk_e_z = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_e_z));
  dk_h_x = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_h_x));
  dk_h_y = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_h_y));
  dk_h_z = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_h_z));

  // initialise the derivative-shift operators for each field component
  init_diff_shift_op(-0.5, dk_e_x, N_e_x);
  init_diff_shift_op(-0.5, dk_e_y, N_e_y);
  init_diff_shift_op(-0.5, dk_e_z, N_e_z);
  init_diff_shift_op(0.5, dk_h_x, N_h_x);
  init_diff_shift_op(0.5, dk_h_y, N_h_y);
  init_diff_shift_op(0.5, dk_h_z, N_h_z);
}

PSTDVariables::~PSTDVariables() {
  if (variables_allocated) {
    fftw_free(dk_e_y);
    fftw_free(dk_e_z);
    fftw_free(dk_h_x);
    fftw_free(dk_h_y);
    fftw_free(dk_h_z);
  }
}
