#include "output_matrices.h"

#include <stdexcept>
#include <algorithm>

#include <spdlog/spdlog.h>

#include "input_output_names.h"
#include "matrix.h"
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

bool OutputMatrices::memory_already_assigned(vector<string> matrix_names) {
  for(string &matrix_name : matrix_names) {
    if (matrix_pointers[index_from_matrix_name(matrix_name)] != nullptr) {
      return true;
    }
  }
  return false;
}
bool OutputMatrices::memory_already_assigned(string matrix_name) {
  if (matrix_pointers[index_from_matrix_name(matrix_name)] != nullptr) {
    return true;
  }
  return false;
}

void OutputMatrices::assign_empty_matrix(std::vector<std::string> matrix_names, mxClassID data_type,
                                         mxComplexity complexity, int ndims) {
  // avoid memory leaks
  error_on_memory_assigned(matrix_names);
  // assign the empty array to all the matrix_names
  vector<int> dims(ndims, 0);
  for (std::string matrix : matrix_names) {
    matrix_pointers[index_from_matrix_name(matrix)] =
            mxCreateNumericArray(ndims, (const mwSize *) dims.data(), data_type, complexity);
  }
}

void OutputMatrices::set_maxresfield(double maxfield, bool overwrite_existing) {
  if (overwrite_existing) {
    if (!memory_already_assigned("maxresfield")) {
      throw runtime_error("Error: maxresfield has no allocated memory, but a write was attempted");
    }
    // overwrite the value currently at maxresfield
  } else {
    // avoid memory leaks
    error_on_memory_assigned("maxresfield");
    // create and allocate new memory
    int dims[2] = {1, 1};
    matrix_pointers[index_from_matrix_name("maxresfield")] =
            mxCreateNumericArray(2, (mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  }
  // assign new maxfield to the output
  *mxGetPr(matrix_pointers[index_from_matrix_name("maxresfield")]) = maxfield;
}

void OutputMatrices::allocate_field_and_labels_memory(bool empty_allocation, int I_tot, int J_tot,
                                                      int K_tot) {
  // matrices we will be assigning in this allocation method
  vector<string> matrices_to_assign = { "Ex_out", "Ey_out", "Ez_out", "Hx_out", "Hy_out", "Hz_out", "x_out", "y_out", "z_out" };

  // create output by reserving memory
  if (empty_allocation) {
    // assign empty matrices if this was requested
    assign_empty_matrix(matrices_to_assign);
  } else {
    // avoid memory leaks
    error_on_memory_assigned(matrices_to_assign);

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

void OutputMatrices::allocate_Id_memory(bool empty_allocation) {
  // create output by reserving memory
  if (empty_allocation) {
    assign_empty_matrix({"Id"});
  } else {
    // avoid memory leaks
    error_on_memory_assigned("Id");

    // create a structure array with two fields
    int dims[2] = {1, 1};
    const char *fieldnames[] = {"Idx", "Idy"};
    matrix_pointers[index_from_matrix_name("Id")] = mxCreateStructArray(2, (const mwSize *) dims, 2, fieldnames);
  }
}

void OutputMatrices::allocate_campssample_memory(int n_vertices, int n_components,
                                                 int n_frequencies) {
  if (n_vertices <= 0) { assign_empty_matrix({"campssample"}, mxDOUBLE_CLASS, mxCOMPLEX, 3); }
  else {
    // avoid memory leaks
    error_on_memory_assigned("campssample");

    // allocate appropriate memory
    int dims[3] = {n_vertices, n_components, n_frequencies};
    matrix_pointers[index_from_matrix_name("campssample")] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  }
}

void OutputMatrices::allocate_interpolation_memory(bool empty_allocation, ElectricField &E, MagneticField &H, Dimension simulation_dimension) {
  // output names that we will be allocating here
  vector<string> interpolated_gridlabels_names = {"x_i", "y_i", "z_i"};
  vector<string> interpolated_E_field_names = {"Ex_i", "Ey_i", "Ez_i"};
  vector<string> interpolated_H_field_names = {"Hx_i", "Hy_i", "Hz_i"};

  if (empty_allocation) {
    // assign empty matrices
    assign_empty_matrix(interpolated_E_field_names);
    assign_empty_matrix(interpolated_H_field_names);
    assign_empty_matrix(interpolated_gridlabels_names);
  } else {
    // avoid memory leaks
    error_on_memory_assigned(interpolated_E_field_names);
    error_on_memory_assigned(interpolated_H_field_names);
    error_on_memory_assigned(interpolated_gridlabels_names);

    // interpolated field memory needs to be assigned
    int i_upper, i_lower, j_upper, j_lower, k_upper, k_lower;
    int outdims[3] = {0, 0, 0};

    // E-field memory assignment
    i_upper = E.I_tot - 2, i_lower = 2;
    j_upper = E.J_tot - 2, j_lower = 2;
    if (simulation_dimension == THREE) {
      k_upper = E.K_tot - 2;
      k_lower = 2;
    } else {
      k_upper = 0;
      k_lower = 0;
    }
    outdims[0] = i_upper - i_lower + 1;
    outdims[1] = j_upper - j_lower + 1;
    if (simulation_dimension == THREE) {
      outdims[2] = k_upper - k_lower + 1;
      outdims[1] = max(1, j_upper - j_lower + 1);
    }
    for(string matrix : interpolated_E_field_names) {
      matrix_pointers[index_from_matrix_name(matrix)] = mxCreateNumericArray(3, (mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
    }
    // H-field memory assignment
    i_upper = H.I_tot - 2, i_lower = 2;
    j_upper = H.J_tot - 2, j_lower = 2;
    if (simulation_dimension == THREE) {
      k_upper = H.K_tot - 2;
      k_lower = 2;
    } else {
      k_upper = 0;
      k_lower = 0;
    }
    outdims[0] = i_upper - i_lower + 1;
    outdims[1] = j_upper - j_lower + 1;
    if (simulation_dimension == THREE) {
      outdims[2] = k_upper - k_lower + 1;
      outdims[1] = max(1, j_upper - j_lower + 1);
    }
    for (string matrix : interpolated_H_field_names) {
      matrix_pointers[index_from_matrix_name(matrix)] = mxCreateNumericArray(3, (mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
    }

    // interpolated gridlabels memory needs to be assigned
    // x-labels are always a fixed size
    int label_dims[2] = {1, E.I_tot - 3};
    matrix_pointers[index_from_matrix_name("x_i")] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);
    // y-labels might recieve nJ < 0, in which case a manual adjustment is made
    label_dims[1] = max(E.J_tot - 3, 1);
    matrix_pointers[index_from_matrix_name("y_i")] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);
    // z-labels depend on the dimensionality of the simulation
    if (simulation_dimension == THREE) {
      label_dims[1] = E.K_tot - 3;
    } else {
      label_dims[1] = 1;
    }
    matrix_pointers[index_from_matrix_name("z_i")] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);
  }
}

void OutputMatrices::allocate_extracted_phasor_memory(bool empty_allocation,
                                                      SurfacePhasors &surface_phasors, mxArray *mx_surface_facets) {
  // output names that we will be assigning to
  vector<string> phasor_matrices = {"vertices", "camplitudes", "facets"};

  if (empty_allocation) {
    assign_empty_matrix(phasor_matrices, mxDOUBLE_CLASS, mxREAL);
  } else {
    // avoid memory leaks
    error_on_memory_assigned(phasor_matrices);

    // pull the memory blocks from the phasors
    matrix_pointers[index_from_matrix_name("vertices")] = surface_phasors.get_vertex_list();
    matrix_pointers[index_from_matrix_name("camplitudes")] =
            surface_phasors.get_mx_surface_amplitudes();
    matrix_pointers[index_from_matrix_name("facets")] = mx_surface_facets;
  }
}

void OutputMatrices::save_outputs(string output_file_name, bool compressed_output) {
  MATFile *output_file = matOpen(output_file_name.c_str(), "w7.3");

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
    // returns 0 if successful and non-zero if an error occured
    int mpv_out = matPutVariable(output_file, matrix_to_write.c_str(), matrix_pointers[index_from_matrix_name(matrix_to_write)]);
    if (mpv_out != 0) {
        // print error information
        FILE *fp = matGetFp(output_file);
        spdlog::info("Could not write array {0:s} to {1:s} ({2:d},{3:d},{4:d})",
                     matrix_to_write.c_str(), output_file_name, mpv_out, feof(fp), ferror(fp));
    }
  }

  // close matfile
  matClose(output_file);
}

OutputMatrices::~OutputMatrices() {
  // search through all possible output matrices
  for(string matrix_name : outputmatrices_all) {
    // check for allocated memory
    if(memory_already_assigned(matrix_name)) {
      // destroy MATLAB structure recursively
      mxDestroyArray(matrix_pointers[index_from_matrix_name(matrix_name)]);
    }
  }
}
