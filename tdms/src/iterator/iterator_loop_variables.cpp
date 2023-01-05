#include "iterator_loop_variables.h"

#include <omp.h>

#include "mesh_base.h"
#include "numerical_derivative.h"

using namespace tdms_phys_constants;
using namespace std;

Iterator_LoopVariables::Iterator_LoopVariables(InputMatrices matrices_from_input_file,
                                               SolverMethod _solver_method,
                                               PreferredInterpolationMethods interpolation_method)
    : Iterator_ObjectsFromInfile(matrices_from_input_file, _solver_method) {

  // set interpolation method for the fields
  set_interpolation_method(interpolation_method);

  // deduce the number of non-pml cells in the z-direction
  K = K_tot - params.pml.Dxl - params.pml.Dxu;

  // deduce refractive index, and print to log
  refind = sqrt(1. / (freespace_Cbx[0] / params.dt * params.delta.dx) / EPSILON0);
  spdlog::info("refind (Refractive index) = {0:e}", refind);

  // setup temporary storgage for detector sensitivity evaluation
  Ex_t = DetectorSensitivityArrays();
  Ey_t = DetectorSensitivityArrays();
  if (params.exdetintegral) {
    int n0 = I_tot - params.pml.Dxl - params.pml.Dxu;
    int n1 = J_tot - params.pml.Dyl - params.pml.Dyu;
    Ex_t.initialise(n1, n0);
    Ey_t.initialise(n1, n0);
  }

  // if running PSTD solver, allocate additional memory for derivative-shift operators now
  if (solver_method == SolverMethod::PseudoSpectral) {
    setup_PSTD_exclusive_variables();
  }

  // initialise the {E,H}_norm variables to an array of zeros
  auto E_norm = (complex<double> *) malloc(f_ex_vec.size() * sizeof(complex<double>));
  auto H_norm = (complex<double> *) malloc(f_ex_vec.size() * sizeof(complex<double>));
  for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
    E_norm[ifx] = 0.;
    H_norm[ifx] = 0.;
  }

  // setup the surface mesh, if it is required
  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    setup_surface_mesh();
  }

  // Now set up the dimensions for the (field) phasor arrays, E and H.
  setup_field_dimensions();
}

void Iterator_LoopVariables::setup_PSTD_exclusive_variables() {
  int max_IJK = E_s.max_IJK_tot(), n_threads = omp_get_max_threads();

  ca_vec.allocate(n_threads, max_IJK + 1);
  cb_vec.allocate(n_threads, max_IJK + 1);
  cc_vec.allocate(n_threads, max_IJK + 1);
  eh_vec.allocate(n_threads, max_IJK + 1);

  // prepare the fourier transform plans
  E_s.initialise_fftw_plan(n_threads, eh_vec);
  H_s.initialise_fftw_plan(n_threads, eh_vec);

  // deduce the number of coefficients in the derivative-shift operator
  N_e_x = I_tot - 1 + 1;
  N_e_y = J_tot - 1 + 1;
  N_e_z = K_tot - 1 + 1;
  N_h_x = I_tot + 1;
  N_h_y = J_tot + 1;
  N_h_z = K_tot + 1;

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

void Iterator_LoopVariables::set_interpolation_method(PreferredInterpolationMethods pim) {
  E_s.set_preferred_interpolation_methods(pim);
  H_s.set_preferred_interpolation_methods(pim);
  E.set_preferred_interpolation_methods(pim);
  H.set_preferred_interpolation_methods(pim);
  E_copy.set_preferred_interpolation_methods(pim);
}

void Iterator_LoopVariables::setup_surface_mesh() {
  // we aren't going to need the facets so create a pointer then delete the memory assigned to it when we're done
  mxArray *temporary_surface_facets;
  if (J_tot == 0) {
    conciseCreateBoundary(cuboid[0], cuboid[1], cuboid[4], cuboid[5], &mx_surface_vertices,
                          &temporary_surface_facets);
  } else {
    conciseTriangulateCuboidSkip(cuboid[0], cuboid[1], cuboid[2], cuboid[3], cuboid[4], cuboid[5],
                                 params.spacing_stride, &mx_surface_vertices,
                                 &temporary_surface_facets);
  }
  // we don't need the facets so destroy the matrix now to save memory
  mxDestroyArray(temporary_surface_facets);

    // now setup the object that handles the surface phasors
  surface_phasors.set_from_matlab_array(mx_surface_vertices, f_ex_vec.size());
  // LEFTOVER comment - unsure of origin and author? (Predates refactor iterator [1/N] (PR #64, 1f4207e25)
  //now need to add a command to update the complex amplitudes
}

void Iterator_LoopVariables::setup_field_dimensions() {
  /*
    Phasors are extracted over the range Dxl + 3 - 1 to I_tot - Dxu - 1 to avoid pml cells - see page III.80 for explanation of the following.

    This has been extended so that interpolation is done at the end of the FDTD run and also to handle the case of when there is no PML in place more appropriately
  */
  E.il = H.il = (params.pml.Dxl) ? params.pml.Dxl + 2 : 0;
  E.iu = H.iu = (params.pml.Dxu) ? I_tot - params.pml.Dxu - 1 : I_tot;
  E.jl = H.jl = (params.pml.Dyl) ? params.pml.Dyl + 2 : 0;
  E.ju = H.ju = (params.pml.Dyu) ? J_tot - params.pml.Dyu - 1 : J_tot;
  E.kl = H.kl = (params.pml.Dzl) ? params.pml.Dzl + 2 : 0;
  E.ku = H.ku = (params.pml.Dzu) ? K_tot - params.pml.Dzu - 1 : K_tot;

  E.I_tot = H.I_tot = E.iu - E.il + 1;
  E.J_tot = H.J_tot = E.ju - E.jl + 1;
  E.K_tot = H.K_tot = E.ku - E.kl + 1;
}

void Iterator_LoopVariables::zero_field_arrays() {
  E.zero();
  H.zero();
  // if converging to steady state, also zero the array that is used to store the phasors from the previous iteration
  if (params.source_mode == SourceMode::steadystate) { E_copy.zero(); }
}

Iterator_LoopVariables::~Iterator_LoopVariables() {
  // if we used the psuedo-spectral method, we need to delete the malloc'd fftw space for the derivative-shift operators
  if (solver_method == SolverMethod::PseudoSpectral) {
    fftw_free(dk_e_x);
    fftw_free(dk_e_y);
    fftw_free(dk_e_z);
    fftw_free(dk_h_x);
    fftw_free(dk_h_y);
    fftw_free(dk_h_z);
  }

  // free {E,H}_norm malloc'd memory
  free(E_norm);
  free(H_norm);

  // free E_copy_data_placeholders if we were running a steady-state simulation
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete) {
    mxDestroyArray(E_copy_data_placeholders[0]);
    mxDestroyArray(E_copy_data_placeholders[1]);
    mxDestroyArray(E_copy_data_placeholders[2]);
  }
}
