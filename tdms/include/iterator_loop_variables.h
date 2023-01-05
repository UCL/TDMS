#pragma once

#include <complex>
#include <fftw3.h>

#include "matrix.h"

#include "arrays.h"
#include "field.h"
#include "globals.h"
#include "input_matrices.h"
#include "iterator_intermediate_matlab_variables.h"
#include "iterator_objects_from_infile.h"
#include "surface_phasors.h"

// Declares all variables used in the main iterator loop that are NOT loaded from the input files
class Iterator_PreLoopDeclarations {
public:
  double refind;//< refractive index of the first layer of the multilayer, or of the bulk of homogeneous

  DetectorSensitivityArrays Ex_t, Ey_t;//< temporary storage for detector sensitivity evaluation

  ElectricField E, E_copy;
  MagneticField H;

  ElectricSplitField E_nm1;
  CurrentDensitySplitField J_c, J_s, J_nm1;

  SurfacePhasors surface_phasors;

  GridLabels output_grid_labels;

  EHVec eh_vec;
  CCoefficientMatrix ca_vec, cb_vec, cc_vec;

  bool is_conductive, is_disp;

  int J_tot_p1_bound, J_tot_bound;

  //refractive index of the first layer of the multilayer, or of the bulk of homogeneous
  double refind;
  double phaseTermE;
  std::complex<double> cphaseTermE;
  double lambda_an_t;

  double Enp1, Jnp1;//< the C and D vars for free space and pml
  double Ca, Cb, Cc;//< used by the interpolation scheme

  double rho;
  double alpha_l, beta_l, gamma_l;
  double kappa_l, sigma_l;

  int i, j, k;
  int n, k_loc;
  int K;//< Number of non-pml cells in the K-direction (K_tot - Dxl - Dxu)

  //these are used for boot strapping. There is currently no way of exporting this.
  double **iwave_lEx_Rbs, **iwave_lEy_Rbs, **iwave_lHx_Rbs, **iwave_lHy_Rbs, **iwave_lEx_Ibs,
          **iwave_lEy_Ibs, **iwave_lHx_Ibs, **iwave_lHy_Ibs;

  std::complex<double> Idxt, Idyt, kprop;
  std::complex<double> *E_norm, *H_norm;
  std::complex<double> **Idx, **Idy;

  /* PSTD storage : variables are not used when using FDTD method */

  // The number of complex fourier coefficients in the derivative-shift operator, for each field component ( N_e_x = # coeffs for Ex, for example)
  int N_e_x, N_e_y, N_e_z, N_h_x, N_h_y, N_h_z;
  // The coefficients of the PSTD derivative shift operator, for each field component ( dk_e_x = coeffs for Ex, for example )
  fftw_complex *dk_e_x, *dk_e_y, *dk_e_z, *dk_h_x, *dk_h_y, *dk_h_z;

  /* MATLAB datatype and pointer declarations */

  mxArray *E_copy_data_placeholders[3];//< Placeholder pointers to the data in the copy of the E-field that checks phasor convergence

  /* MATLAB datatype pointers that are directly assigned to outputs, so are not free'd on destruction */

  mxArray *mx_surface_vertices = nullptr;//< Vertices on the user-defined surface
};

// Handles all the variables declared before the main iterator loop, and their initalisation if appropriate
class Iterator_LoopVariables : public Iterator_ObjectsFromInfile, Iterator_PreLoopDeclarations {
private:
  /**
   * @brief  Set the interpolation method for each of the fields we are interpolating
   *
   * @param pim The interpolation method to use
   */
  void set_interpolation_method(PreferredInterpolationMethods pim);
  // Setup the PSTD derivative-shift operator and related variables
  void setup_PSTD_exclusive_variables();
  // Sets up the mesh for the user-specified surface, over which to extract phasors
  void setup_surface_mesh();
  // Set the dimensions of the field (phasor) arrays for this simulation
  void setup_field_dimensions();
  // Zero the field arrays
  void zero_field_arrays();

public:
    Iterator_LoopVariables(InputMatrices matrices_from_input_file, SolverMethod _solver_method, PreferredInterpolationMethods interpolation_method);
  /**
   * @brief Sets the pointers for the output field and grid_labels data
   *
   * @param output_pointers The array of output pointers to assign to
   */
  void link_fields_and_labels_to_output_pointers(mxArray *output_pointers[]);
  /**
   * @brief Sets the pointers for the output Id object and allocates the appropriate class members
   *
   * @param output_pointers The array of output pointers to assign to
   */
  void link_id_to_output_pointers(mxArray *output_pointers[], Iterator_IntermediateMATLABVariables &iMVars);

    ~Iterator_LoopVariables();
};
