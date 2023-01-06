#pragma once

#include <complex>
#include <fftw3.h>

#include "matrix.h"

#include "arrays.h"
#include "field.h"
#include "globals.h"
#include "input_matrices.h"
#include "iterator_objects_from_infile.h"
#include "surface_phasors.h"

/* Declares all variables used in the main iterator loop that are NOT loaded from the input files
This class is NOT intended to see any standalone use - it exists soley to be inherited by Iterator_LoopVariables.

The default destructor for this class can handle cleaning up all of our custom objects.
No members of this class have memory malloc'd, so there is no need for manual memory freeing in a destructor.
 */
class Iterator_PreLoopDeclarations {
public:
  DetectorSensitivityArrays Ex_t, Ey_t;//< temporary storage for detector sensitivity evaluation

  ElectricField E, E_copy;
  MagneticField H;

  ElectricSplitField E_nm1;
  CurrentDensitySplitField J_c;//< The per-cell ( current density or conductivity ? ) of the material
  CurrentDensitySplitField J_s, J_nm1;

  SurfacePhasors surface_phasors;

  GridLabels output_grid_labels;

  EHVec eh_vec;
  CCoefficientMatrix ca_vec, cb_vec, cc_vec;

  bool is_conductive, is_disp;

  int J_tot_p1_bound, J_tot_bound;

  double refind;//< refractive index of the first layer of the multilayer, or of the bulk of homogeneous
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

  std::complex<double> Idxt, Idyt, kprop;
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
  // Determine the dispersive properties and setup the corresponding arrays
  void setup_dispersive_properties();
  /**
   * @brief Determine if we have a dispersive medium by searching for non-zero attenuation constants in gamma.
   *
   * @param non_zero_tol Tolerance for an attenuation constant being "non-zero"
   * @return true, we have a dispersive medium
   * @return false, we do not have a dispersive medium
   */
  bool is_dispersive(double non_zero_tol = 1e-15);

  // Zero the field arrays
  void zero_field_arrays();
  /**
   * @brief Zero all elements of a n_dim1-by-n_dim2 array that has been cast to a MATLAB array.
   *
   * The array is assumed to be indexed by cast_array[j][i], where
   * 0 <= j < n_dim2,
   * 0 <= i < n_dim1.
   *
   * @param cast_array The array whose elements are to be set to zero
   * @param n_dim1,n_dim2 The number of elements in each axis
   */
  void zero_cast_array(double **cast_array, int n_dim1, int n_dim2);

protected:
  /* INTERMEDIATE MATLAB VARIABLES/POINTERS */
  /* MATLAB DANGER ZONE: Some pointers are created via malloc-like commands, but are linked to would-be outputs. These are destroyed when their parent arrays are destroyed via mxFree, however this can't be done by this class since we need to pass back the output data before we delete it. Need to check that mxFree is called on plhs after writing outputs! */

  mxArray *E_copy_data_placeholders
          [3];//< Placeholder pointers to the data in the copy of the E-field that checks phasor convergence

  mxArray *mx_surface_vertices = nullptr;//< Vertices on the user-defined surface

  mxArray *mx_Idx, *mx_Idy;//< Holds the arrays in the Idx and Idy fields of plhs[26] (Id output) MIGHT BE ABLE TO GET AWAY WITH NOT DECLARING THESE, AND JUST DIRECTLY ASSIGNING TO PLHS[26]. INVESTIGATE AFTER REFACTOR

  /* OTHER VARIABLES THAT NEED TO BE FREE'D */

  std::complex<double> *E_norm, *H_norm;
  std::complex<double> **Idx, **Idy;

  double **Idx_re, **Idy_re, **Idx_im,
          **Idy_im;//< Point to the real (re) and imaginary (im) parts of the data in the Id output

  // Frees any non-output-required memory assigned to the Id variables
  void free_Id_memory();

  //these are used for fdtd boot strapping. There is currently no way of exporting this.
  double **iwave_lEx_Rbs, **iwave_lEy_Rbs, **iwave_lHx_Rbs, **iwave_lHy_Rbs, **iwave_lEx_Ibs,
          **iwave_lEy_Ibs, **iwave_lHx_Ibs, **iwave_lHy_Ibs;

  // Frees any non-output-required memory assigned to the iwave_l variables used in bootstrapping
  void free_iwave_memory();

public:
  /* PSTD EXCLUSIVE VARIABLES (REQUIRES MEMORY MANAGEMENT)
  These variables are not used when using FDTD method */

  // The number of complex fourier coefficients in the derivative-shift operator, for each field component ( N_e_x = # coeffs for Ex, for example)
  int N_e_x, N_e_y, N_e_z, N_h_x, N_h_y, N_h_z;
  // The coefficients of the PSTD derivative shift operator, for each field component ( dk_e_x = coeffs for Ex, for example )
  fftw_complex *dk_e_x, *dk_e_y, *dk_e_z, *dk_h_x, *dk_h_y, *dk_h_z;

  Iterator_LoopVariables(InputMatrices matrices_from_input_file, SolverMethod _solver_method, PreferredInterpolationMethods interpolation_method);

  /**
   * @brief Sets the pointers for the output field and grid_labels data
   *
   * @param output_pointers The array of output pointers to assign to
   */
  void link_fields_and_labels(mxArray *output_pointers[]);

  /**
   * @brief Sets the pointers for the output Id object and allocates the appropriate class members
   *
   * @param output_pointers The array of output pointers to assign to
   */
  void link_id(mxArray *output_pointers[]);

  /**
   * @brief Setup the phasor arrays for storing the fdtd version of the input fields.
   *
   * These will be used in a boot strapping procedure. Calculated over a complete
    xy-plane.
   *
   * @param output_pointers The array of output pointers to assign to
   */
  void link_fdtd_phasor_arrays(mxArray *output_pointers[]);

  ~Iterator_LoopVariables();
};
