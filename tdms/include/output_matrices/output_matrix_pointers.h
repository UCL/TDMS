/**
 * @file output_matrix_pointers.h
 * @brief Declares a class that allows for memory-management of pointers to the output arrays (OutputMatrixPointers)
 */
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "input_output_names.h"
#include "matrix.h"

/**
 * @brief Provides human-readable and safe management of the pointers to the output arrays, and the associated memory.
 *
 * Instances can be referred using instance[matrix_name] to obtain the pointer to the output matrix matrix_name. However, memory management checks exist to prevent us from leaking memory through blase assignment.
 */
class OutputMatrixPointers {
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

public:
  /**
   * @brief Fetch a (pointer to a) matrix by name.
   *
   * @param matrix_name Name of the matrix to fetch the pointer to
   * @return mxArray* Pointer to the corresponding MATLAB array
   */
  mxArray *&operator[](const std::string &matrix_name) {
    return matrix_pointers[index_from_matrix_name(matrix_name)];
  }

  /**
   * @brief Return true if the matrix_pointer corresponding to any one of the matrix names passed already points to allocated memory.
   *
   * @param matrix_names List of output matrices to check.
   */
  bool memory_already_assigned(const std::vector<std::string> &matrix_names);

  /**
   * @brief Throw an error if the matrix_pointer corresponding to any one of the matrix names passed already points to allocated memory.
   *
   * This function exists to avoid memory leaks - it is possible for us to reassign the matrix_pointers after mxCreate-ing arrays. To avoid this, we always validate we have not done this before attempting to assign memory.
   *
   * @param matrix_names List of output matrices to check.
   */
  void error_on_memory_assigned(const std::vector<std::string> &matrix_names) {
    if (memory_already_assigned(matrix_names)) {
      throw std::runtime_error("Reassigning output pointer will cause memory leak - aborting.");
    }
  }

  /**
   * @brief Assign the matrices provided as empty arrays.
   *
   * Outputs are typically set to empty arrays when the user has requested the simulation be run in such a way that the output is not generated, or not requested. We still check for memory leaks before assignment as creating an empty array can still leave MATLAB memory not-cleaned-up.
   *
   * @param matrix_names The list of matrix names to set to be empty.
   * @param data_type The data type to assign to the empty array (needed for compatibility issues)
   * @param is_complex Whether the data is complex (mxCOMPLEX) or real (mxREAL)
   * @param ndims Number of dimensions in the empty array to assign
   */
  void assign_empty_matrix(const std::vector<std::string> &matrix_names,
                           mxClassID data_type = mxDOUBLE_CLASS,
                           mxComplexity complexity = mxCOMPLEX, int ndims = 2);

  ~OutputMatrixPointers();
};
