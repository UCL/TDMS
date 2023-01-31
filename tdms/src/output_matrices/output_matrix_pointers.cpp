#include "output_matrices/output_matrix_pointers.h"

#include <algorithm>

using namespace tdms_matrix_names;
using namespace std;

int OutputMatrixPointers::index_from_matrix_name(const string &matrix_name) {
  auto position = find(outputmatrices_all.begin(), outputmatrices_all.end(), matrix_name);
  if (position == outputmatrices_all.end()) {
    // could not find the matrix name in the list of expected input matrices
    throw runtime_error(matrix_name + " not found in outputmatrices_all");
  }
  return distance(outputmatrices_all.begin(), position);
}

bool OutputMatrixPointers::memory_already_assigned(const vector<string> &matrix_names) {
  for (const string &matrix_name : matrix_names) {
    if (matrix_pointers[index_from_matrix_name(matrix_name)] != nullptr) { return true; }
  }
  return false;
}

void OutputMatrixPointers::assign_empty_matrix(const vector<string> &matrix_names,
                                               mxClassID data_type, mxComplexity complexity,
                                               int ndims) {
  // avoid memory leaks
  error_on_memory_assigned(matrix_names);
  // assign the empty array to all the matrix_names
  vector<int> dims(ndims, 0);
  for (string matrix : matrix_names) {
    matrix_pointers[index_from_matrix_name(matrix)] =
            mxCreateNumericArray(ndims, (const mwSize *) dims.data(), data_type, complexity);
  }
}

OutputMatrixPointers::~OutputMatrixPointers() {
  // search through all possible output matrices
  for (string matrix_name : outputmatrices_all) {
    // check for allocated memory
    if (memory_already_assigned({matrix_name})) {
      // destroy MATLAB structure recursively
      mxDestroyArray(matrix_pointers[index_from_matrix_name(matrix_name)]);
    }
  }
}
