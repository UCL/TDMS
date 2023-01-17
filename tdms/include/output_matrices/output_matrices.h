/**
 * @file output_matrices.h
 * @brief Class handling the data structures to be written to the output file.
 */
#pragma once

#include <string>

#include "cell_coordinate.h"
#include "field.h"
#include "fieldsample.h"
#include "grid_labels.h"
#include "id_variables.h"
#include "matrix.h"
#include "output_matrix_pointers.h"
#include "simulation_parameters.h"
#include "surface_phasors.h"
#include "vertex_phasors.h"

/**
 * @brief Handles the data that is written to the output file, and the associated native C++ classes and datatypes.
 *
 * Some of the C++ datatypes are updated in the main iterative loop, whilst others are obtained using the data from the main loop. This class handles both cases, and uses pointers to track where the MATLAB-structures to be written out are stored.
 *
 * Upon writing the output, an instance of this class can go out of scope and will take with it any reserved memory for MATLAB structures and native C++ arrays of dynamic size.
 */
class OutputMatrices {
private:
  // Pointers to arrays in C++ that will be populated by pointers to the output data
  OutputMatrixPointers output_arrays;

  IJKDims n_Yee_cells;//< Number of Yee cells in each axial direction for the simulation (dictated by E_s split-field)

  /**
   * @brief Computes the field values at the centre of the Yee cells, and the corresponding spatial coordinates of these locations.
   *
   * @param dimension Dimensionality of the simulation: THREE, TE, or TM.
   */
  void compute_interpolated_fields(Dimension dimension);

public:
  OutputMatrices() = default;

  /**
   * @brief Record the number of Yee cells in each axial direction from the split-electric field that was input.
   *
   * @param IJK_tot Number of Yee cells in the {I,J,K} directions
   */
  void set_n_Yee_cells(IJKDims IJK_tot) { n_Yee_cells = IJK_tot; }

  /**
   * @brief Fetch a (pointer to a) MATLAB output matrix by name.
   *
   * @param matrix_name Name of the matrix to fetch the pointer to
   * @return mxArray* Pointer to the corresponding MATLAB array
   */
  mxArray *&operator[](const std::string &matrix_name) {
    return output_arrays[matrix_name];
  }

  IDVariables ID;//< The ID output, and associated variables

  /**
   * @brief Create MATLAB memory for the Id structure array.
   *
   * @param empty_allocation If true, empty arrays will be allocated
   * @param n_frequencies Number of frequencies that we're extracting at
   * @param n_det_modes D_tilde.num_det_modes()
   */
  void setup_Id(bool empty_allocation, int n_frequencies, int n_det_modes);

  ElectricField E;//< Electric field and phasors at the Yee cell positions
  MagneticField H;//< Magnetic field and phasors at the Yee cell positions
  GridLabels output_grid_labels;//< Co-ordinates (spatial positions) of the field values

  /**
   * @brief Set the up E, H, and output_grid_labels outputs
   *
   * @param params The simulation parameters for this run
   * @param input_grid_labels The grid labels obtained from the input file
   */
  void setup_EH_and_gridlabels(SimulationParameters params, GridLabels input_grid_labels);
  /**
   * @brief Get the dimensions of the electric (and magnetic) field.
   *
   * @return IJKDims The dimensions of the electric (and magnetic) field
   */
  IJKDims get_E_dimensions() {
    return IJKDims(E.I_tot, E.J_tot, E.K_tot);
  }

  SurfacePhasors surface_phasors;//< Phasors extracted over the user-specified surface

  /**
   * @brief Create MATLAB memory for the surface phasor outputs
   *
   * @param empty_allocation If true, empty arrays will be allocated
   * @param mx_surface_facets The surface facets to write out
   */
  void assign_surface_phasor_outputs(bool empty_allocation, mxArray *mx_surface_facets);

  VertexPhasors vertex_phasors;//< Phasors extracted at user-specified vertices

  /**
   * @brief Setup the object that will handle extraction of the phasors at the vertices
   *
   * @param vp_ptr Pointer to the input array containing vertex data
   * @param n_frequencies The number of frequencies we are extracting at
   */
  void setup_vertex_phasors(const mxArray *vp_ptr, int n_frequencies);

  FieldSample fieldsample;//< E,H split-field values sampled at user-specified locations and frequencies

  /**
   * @brief Setup the object that handles extraction of the field values at the sample positions
   *
   * @param fieldsample_input_data Pointer to the input array containing sample positions
   */
  void setup_fieldsample(const mxArray *fieldsample_input_data);

  /**
   * @brief Set the maxresfield object. If memory has not been reserved yet, it can be assigned here.
   *
   * @param maxfield Value to write to output memory
   * @param overwrite_existing If true, overwrite existing value (do not assign new memory). Otherwise, create new memory for the maxresfield output.
   */
  void set_maxresfield(double maxfield, bool overwrite_existing);

  GridLabels interp_output_grid_labels;//< Holds the spatial co-ordinates of the interpolated fields

  /**
   * @brief Create MATLAB memory for the gridlabels for the interpolated fields
   *
   * @param empty_allocation If true, empty arrays will be allocated
   * @param E,H The {electric,magnetic} field that will be interpolated
   * @param simulation_dimension Whether we are running a 3D, TE, or TM simulation
   */
  void setup_interpolation_outputs(SimulationParameters params);

  /**
   * @brief Save the output matrices to the output file
   *
   * @param output_file_name The file to write the simulation outputs to
   * @param compressed_output If true, write compressed output (do not write facets and vertices)
   */
  void save_outputs(std::string output_file_name, bool compressed_output = false);
};
