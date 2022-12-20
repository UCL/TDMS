#pragma once

#include "mat_io.h"
#include "input_output_names.h"

int index_from_matrix_name(const char *matrix_name);

class InputMatrices {
private:
  // Pointers to arrays in C++ that will be populated by the MATLAB matrices
  const mxArray *matrix_pointers[NMATRICES];

  void assign_matrix_pointers(MatrixCollection &expected, MatFileMatrixCollection &actual);

public:
    InputMatrices() = default;

    const mxArray *operator[](int index) {return matrix_pointers[index]; };
    const mxArray *operator[](const char *matrix_name) { return matrix_pointers[index_from_matrix_name(matrix_name)]; };

    void set_matrix_pointer(int index, mxArray *new_ptr) { matrix_pointers[index] = new_ptr; };
    void set_matrix_pointer(const char *matrix_name, mxArray *new_ptr) {
      set_matrix_pointer(index_from_matrix_name(matrix_name), new_ptr);
    }

    void set_from_input_file(const char *mat_filename, const char *gridfile="");
};
