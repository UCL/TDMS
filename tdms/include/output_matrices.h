/**
 * @file output_matrices.h
 * @brief Class handling the data structures to be written to the output file.
 */
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "input_output_names.h"
#include "field.h"
#include "matrix.h"
#include "matrix_collection.h"
#include "simulation_parameters.h"
#include "surface_phasors.h"

class OutputMatrices {
private:
  // Pointers to arrays in C++ that will be populated by pointers to the output data
  mxArray *matrix_pointers[NOUTMATRICES_WRITE_ALL] = {nullptr};

  /**
   * @brief Fetches the index of matrix_name in the tdms_matrix_names::matrixnames array
   *
   * @param matrix_name The matrix name to locate
   * @return int The corresponding index
   */
  int index_from_matrix_name(const std::string &matrix_name);

  /**
   * @brief Return true if the matrix_pointer corresponding to any one of the matrix names passed already points to allocated memory.
   *
   * @param matrix_names List of output matrices to check.
   */
  bool memory_already_assigned(std::vector<std::string> matrix_names);
  /*! @copydoc memory_already_assigned */
  bool memory_already_assigned(std::string matrix_name);
  /**
   * @brief Throw an error if the matrix_pointer corresponding to any one of the matrix names passed already points to allocated memory.
   *
   * This function exists to avoid memory leaks - it is possible for us to reassign the matrix_pointers after mxCreate-ing arrays. To avoid this, we always validate we have not done this before attempting to assign memory.
   *
   * @param matrix_names List of output matrices to check.
   */
  void error_on_memory_assigned(std::vector<std::string> matrix_names) {
    if (memory_already_assigned(matrix_names)) {
      throw std::runtime_error("Reassigning output pointer will cause memory leak - aborting.");
    }
  }
  /*! @copydoc error_on_memory_assigned */
  void error_on_memory_assigned(std::string matrix_name) {
    if (memory_already_assigned(matrix_name)) {
      throw std::runtime_error("Reassigning output pointer will cause memory leak - aborting.");
    }
  }
  /**
   * @brief Assign the output matrices provided as empty arrays.
   *
   * Outputs are typically set to empty arrays when the user has requested the simulation be run in such a way that the output is not generated, or not requested. We still check for memory leaks before assignment as creating an empty array can still leave MATLAB memory not-cleaned-up.
   *
   * @param matrix_names The list of matrix names to set to be empty.
   * @param data_type The data type to assign to the empty array (needed for compatability issues)
   * @param is_complex Whether the data is complex (mxCOMPLEX) or real (mxREAL)
   * @param ndims Number of dimensions in the empty array to assign
   */
  void assign_empty_matrix(std::vector<std::string> matrix_names, mxClassID data_type = mxDOUBLE_CLASS, mxComplexity complexity = mxCOMPLEX, int ndims = 2);

public:
  OutputMatrices() = default;

  /**
   * @brief Fetch a (pointer to a) MATLAB output matrix by index reference.
   *
   * Index-references can be computed through index_from_matrix_name. They are ordered in the same way as the names in tdms_matrix_names::matrixnames.
   *
   * @param index The index to fetch
   * @return mxArray* Pointer to the corresponding MATLAB array
   */
  mxArray* &operator[](int index) { return matrix_pointers[index]; }
  /**
   * @brief Fetch a (pointer to a) MATLAB output matrix by name.
   *
   * @param matrix_name Name of the matrix to fetch the pointer to
   * @return mxArray* Pointer to the corresponding MATLAB array
   */
  mxArray* &operator[](const std::string &matrix_name) {
    return matrix_pointers[index_from_matrix_name(matrix_name)];
  }

  /**
   * @brief Set the matrix pointer with the given index.
   *
   * Not recommended unless you're certain you are setting the correct array pointer! Use the overload which takes the array name as an input if you are uncertain you are setting the correct pointer.
   *
   * This function mainly sees use during destruction (and output writing) when we just want to iterate over all the matrices that we are dealing with.
   *
   * @param index The index in matrix_pointers to set
   * @param new_ptr The address we want to assign
   */
  void set_matrix_pointer(int index, mxArray *new_ptr) { matrix_pointers[index] = new_ptr; }
  /**
   * @brief Set the pointer to the named matrix
   *
   * @param matrix_name The name of the matrix that we want to set the pointer to
   * @param new_ptr The address we want to assign
   */
  void set_matrix_pointer(const std::string &matrix_name, mxArray *new_ptr) {
    set_matrix_pointer(index_from_matrix_name(matrix_name), new_ptr);
  }
  /**
   * @brief Set the maxresfield object. If memory has not been reserved yet, it can be assigned here.
   *
   * @param maxfield Value to write to output memory
   * @param overwrite_existing If true, overwrite existing value (do not assign new memory). Otherwise, create new memory for the maxresfield output.
   */
  void set_maxresfield(double maxfield, bool overwrite_existing);

  /**
   * @brief Create MATLAB memory for the field and gridlabel outputs.
   *
   * @param empty_allocation If true, empty arrays will be allocated
   * @param I_tot,J_tot,K_tot The dimensions of the field components (corresponding gridlabel dimensions will be inferred)
   */
  void allocate_field_and_labels_memory(bool empty_allocation, int I_tot = 0, int J_tot = 0, int K_tot = 0);
  /**
   * @brief Create MATLAB memory for the Id structure array.
   *
   * @param empty_allocation If true, empty arrays will be allocated
   */
  void allocate_Id_memory(bool empty_allocation);
  /**
   * @brief Create MATLAB memory for the campssample array.
   *
   * @param n_vertices Number of vertices we are extracting at. If this is 0, then allocate the empty array
   * @param n_components Number of amplitude components we are extracting
   * @param n_frequencies Number of frequencies we are extracting at
   */
  void allocate_campssample_memory(int n_vertices, int n_components, int n_frequencies);
  /**
   * @brief Create MATLAB memory for the gridlabels for the interpolated fields
   *
   * @param empty_allocation If true, empty arrays will be allocated
   * @param E,H The {electric,magnetic} field that will be interpolated
   * @param simulation_dimension Whether we are running a 3D, TE, or TM simulation
   */
  void allocate_interpolation_memory(bool empty_allocation, ElectricField &E, MagneticField &H, Dimension simulation_dimension = THREE);
  /**
   * @brief Create MATLAB memory for the (surface) phasor outputs
   *
   * @param empty_allocation If true, empty arrays will be allocated
   * @param surface_phasors The surface phasor object handling phasor extraction
   * @param mx_surface_facets The surface facets to write out
   */
  void allocate_extracted_phasor_memory(bool empty_allocation, SurfacePhasors &surface_phasors, mxArray *mx_surface_facets);

  /**
   * @brief Save the output matrices to the output file
   *
   * @param output_file_name The file to write the simulation outputs to
   * @param compressed_output If true, write compressed output (do not write facets and vertices)
   */
  void save_outputs(std::string output_file_name, bool compressed_output = false);

  ~OutputMatrices();
};
