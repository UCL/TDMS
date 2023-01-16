#pragma once

#include "arrays.h"
#include "field.h"
#include "grid_labels.h"
#include "surface_phasors.h"

class LoopVariables {
private:
  /**
   * @brief Set the interpolation method for each of the fields we are interpolating
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

public:
  DetectorSensitivityArrays Ex_t, Ey_t;//< temporary storage for detector sensitivity evaluation

  ElectricField E, E_copy;
  MagneticField H;

  ElectricSplitField E_nm1;
  CurrentDensitySplitField
          J_c;//< The per-cell ( current density or conductivity ? ) of the material
  CurrentDensitySplitField J_s, J_nm1;

  SurfacePhasors surface_phasors;

  GridLabels output_grid_labels;

  EHVec eh_vec;
  CCoefficientMatrix ca_vec, cb_vec, cc_vec;

  bool is_conductive, is_disp;

  int J_tot_p1_bound, J_tot_bound;

  double refind;//< refractive index of the first layer of the multilayer, or of the bulk of homogeneous

  int K;//< Number of non-pml cells in the K-direction (K_tot - Dxl - Dxu)
};
