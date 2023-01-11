#include "output_matrices.h"

#include <stdexcept>
#include <algorithm>

#include <spdlog/spdlog.h>

#include "input_output_names.h"
#include "utils.h"

using namespace tdms_matrix_names;
using namespace std;

int OutputMatrices::index_from_matrix_name(const std::string &matrix_name) {
  auto position = find(outputmatrices_all.begin(), outputmatrices_all.end(), matrix_name);
  if (position == outputmatrices_all.end()) {
    // could not find the matrix name in the list of expected input matrices
    throw runtime_error(matrix_name + " not found in outputmatrices_all");
  }
  return distance(outputmatrices_all.begin(), position);
}

bool OutputMatrices::memory_already_assigned(std::vector<std::string> matrix_names) {
  for(string &matrix_name : matrix_names) {
    if (matrix_pointers[index_from_matrix_name(matrix_name)] != nullptr) {
      spdlog::info(matrix_name + " is already assigned!");
      return true;
    }
  }
  return false;
}

void OutputMatrices::assign_empty_matrix(std::vector<std::string> matrix_names, mxClassID data_type,
                                         mxComplexity complexity) {
  // avoid memory leaks
  error_on_memory_assigned(matrix_names);
  // assign the empty array to all the matrix_names
  int dims[2] = {0, 0};
  for (std::string matrix : matrix_names) {
    matrix_pointers[index_from_matrix_name(matrix)] =
            mxCreateNumericArray(2, (const mwSize *) dims, data_type, complexity);
  }
}

void OutputMatrices::allocate_field_and_labels_memory(bool empty_allocation, int I_tot, int J_tot, int K_tot) {
  // matrices we will be assigning in this allocation method
  vector<string> matrices_to_assign = { "Ex_out", "Ey_out", "Ez_out", "Hx_out", "Hy_out", "Hz_out", "x_out", "y_out", "z_out"};
  // avoid memory leaks
  error_on_memory_assigned(matrices_to_assign);

  // create output by reserving memory
  if (empty_allocation) {
    // assign empty matrices if this was requested
    assign_empty_matrix(matrices_to_assign);
  } else {
    // initialise to actual, proper arrays
    int dims[3] = {I_tot, J_tot, K_tot};
    spdlog::info("dims: ({0:d},{1:d},{2:d})", dims[0], dims[1], dims[2]);

    // create MATLAB data storage for the field-component outputs
    matrix_pointers[index_from_matrix_name("Ex_out")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    matrix_pointers[index_from_matrix_name("Ey_out")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    matrix_pointers[index_from_matrix_name("Ez_out")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    matrix_pointers[index_from_matrix_name("Hx_out")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    matrix_pointers[index_from_matrix_name("Hy_out")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    matrix_pointers[index_from_matrix_name("Hz_out")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    // create MATLAB data storage for the gridlabel outputs
    int label_dims_x[2] = {1, I_tot};
    int label_dims_y[2] = {1, J_tot};
    int label_dims_z[2] = {1, K_tot};
    matrix_pointers[index_from_matrix_name("x_out")] = mxCreateNumericArray(2, (const mwSize *) label_dims_x, mxDOUBLE_CLASS, mxREAL);
    matrix_pointers[index_from_matrix_name("y_out")] = mxCreateNumericArray(2, (const mwSize *) label_dims_y, mxDOUBLE_CLASS, mxREAL);
    matrix_pointers[index_from_matrix_name("z_out")] = mxCreateNumericArray(2, (const mwSize *) label_dims_z, mxDOUBLE_CLASS, mxREAL);
  }
}

void OutputMatrices::save_outputs(string output_file_name, bool compressed_output) {
  auto output_file = matOpen(output_file_name.c_str(), "w7.3");

  // check output file was opened successfully
  if (output_file == nullptr) {
    throw runtime_error("Unable to open output file" + string(output_file_name));
  }

  // if we don't want a compressed output we save additional matrices
  vector<string> output_matrices_requested;
  if (compressed_output) {
    output_matrices_requested = outputmatrices;
  } else {
    output_matrices_requested = outputmatrices_all;
  }
  // iterate through the matrices we want to save, setting names and placing them into the matfile
  for(string matrix_to_write : output_matrices_requested) {
    auto mpv_out = matPutVariable(output_file, matrix_to_write.c_str(), (mxArray *) operator[](matrix_to_write));
    if (mpv_out) {
        auto fp = matGetFp(output_file);
        spdlog::info("Could not write array {0:s} to {1:s} ({2:d},{3:d},{4:d})",
                     matrix_to_write.c_str(), output_file_name, mpv_out, feof(fp), ferror(fp));
    }
  }

  // close matfile
  matClose(output_file);
}
