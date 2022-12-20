/**
 * @file input_matrices.h
 * @brief Class handling the input matrices that are read in from the input file
 */
#pragma once

#include "input_output_names.h"
#include "mat_io.h"
#include "matrix_collection.h"

/**
 * @brief Fetches the index of matrix_name in the tdms_matrix_names::matrixnames array
 *
 * @param matrix_name The matrix name to locate
 * @return int The corresponding index
 */
int index_from_matrix_name(const char *matrix_name);

class InputMatrices {
private:
  // Pointers to arrays in C++ that will be populated by the MATLAB matrices (default to nullptrs)
  const mxArray *matrix_pointers[NMATRICES];

  /**
   * @brief Assigns pointers to the matrices in an input file, based on those we are expecting to recieve.
   *
   * MatrixCollection is a list of names of matrices we expect to be in the MatFileMatrixCollection, based on whether we are reading from a gridfile, input_file without a gridfile, or an input file that came with a gridfile.
   *
   * @param expected The matrices we expect to be able to pull from actual
   * @param actual The matrices saved in an input file, to be loaded
   */
  void assign_matrix_pointers(MatrixCollection &expected, MatFileMatrixCollection &actual);

public:
  InputMatrices() = default;

  /**
     * @brief Fetch a (pointer to a) MATLAB input matrix by index reference.
     *
     * Index-references can be computed through index_from_matrix_name. They are ordered in the same way as the names in tdms_matrix_names::matrixnames.
     *
     * @param index The index to fetch
     * @return const mxArray* Pointer to the corresponding MATLAB array
     */
  const mxArray *operator[](int index) { return matrix_pointers[index]; }
  /**
     * @brief Fetch a (pointer to a) MATLAB input matrix by name.
     *
     * @param matrix_name Name of the matrix to fetch the pointer to
     * @return const mxArray* Pointer to the corresponding MATLAB array
     */
  const mxArray *operator[](const char *matrix_name) {
    return matrix_pointers[index_from_matrix_name(matrix_name)];
  }

  /**
     * @brief Set the pointer with the given index
     */
  void set_matrix_pointer(int index, mxArray *new_ptr) { matrix_pointers[index] = new_ptr; }
  /**
     * @brief Set the pointer to the named matrix
     */
  void set_matrix_pointer(const char *matrix_name, mxArray *new_ptr) {
    set_matrix_pointer(index_from_matrix_name(matrix_name), new_ptr);
  }

  /**
   * @brief Open the input mat file, load the matrices, and setup pointers to the matrices
   *
   * @param mat_filename The MATLAB filename
   */
  void set_from_input_file(const char *mat_filename);
  /**
   * @brief Open the input mat file, load the matrices, and setup pointers to the matrices
   *
   * @param mat_filename The MATLAB filename
   * @param gridfile The additional gridfile
   */
  void set_from_input_file(const char *mat_filename, const char *gridfile);
};
