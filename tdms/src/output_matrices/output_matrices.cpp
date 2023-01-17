#include "output_matrices.h"

#include <spdlog/spdlog.h>

using namespace tdms_matrix_names;
using namespace std;

void OutputMatrices::compute_interpolated_fields(Dimension dimension) {
  // link the interpolated gridlabels to the respective outputs
  interp_output_grid_labels.x = mxGetPr(output_arrays["x_i"]);
  interp_output_grid_labels.y = mxGetPr(output_arrays["y_i"]);
  interp_output_grid_labels.z = mxGetPr(output_arrays["z_i"]);

  // depending on whether we are in 3D or not, the Yee cell range we interpolate over changes, which we address here
  int i_upper = E.I_tot - 2, i_lower = 2;
  int j_upper = E.J_tot - 2, j_lower = 2;
  int k_upper = E.K_tot - 2, k_lower = 2;
  if (dimension != THREE) {
    // change K-range in a non-3D simulation
    k_upper = k_lower = 0;
  }

  // now interpolate over the extracted phasors, and set the interpolated gridlabels
  E.interpolate_over_range(output_arrays["Ex_i"], output_arrays["Ey_i"], output_arrays["Ez_i"],
                           i_lower, i_upper, j_lower, j_upper, k_lower, k_upper, dimension);
  H.interpolate_over_range(output_arrays["Hx_i"], output_arrays["Hy_i"], output_arrays["Hz_i"],
                           i_lower, i_upper, j_lower, j_upper, k_lower, k_upper, dimension);
  interp_output_grid_labels.initialise_from(output_grid_labels, i_lower, i_upper, j_lower, j_upper,
                                            k_lower, k_upper);
}

void OutputMatrices::setup_vertex_phasors(const mxArray *vp_ptr, int n_frequencies) {
  // setup from the input provided
  vertex_phasors.set_from(vp_ptr);

  // avoid memory leaks
  output_arrays.error_on_memory_assigned({"campssample"});

  // assign memory if we have been asked to extract phasors at specific vertices
  if (vertex_phasors.there_are_vertices_to_extract_at()) {
    vertex_phasors.setup_complex_amplitude_arrays(n_frequencies);
    output_arrays["campssample"] = vertex_phasors.get_mx_camplitudes();
  } else {
    output_arrays.assign_empty_matrix({"campssample"}, mxDOUBLE_CLASS, mxCOMPLEX, 3);
  }
}

void OutputMatrices::setup_fieldsample(const mxArray *fieldsample_input_data) {
  fieldsample.set_from(fieldsample_input_data);
  output_arrays["fieldsample"] = fieldsample.mx;
}

void OutputMatrices::set_maxresfield(double maxfield, bool overwrite_existing) {
  if (overwrite_existing) {
    if (!output_arrays.memory_already_assigned({"maxresfield"})) {
      throw runtime_error("Error: maxresfield has no allocated memory, but a write was attempted");
    }
    // overwrite the value currently at maxresfield
  } else {
    // avoid memory leaks
    output_arrays.error_on_memory_assigned({"maxresfield"});
    // create and allocate new memory
    int dims[2] = {1, 1};
    output_arrays["maxresfield"] =
            mxCreateNumericArray(2, (mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  }
  // assign new maxfield to the output
  *mxGetPr(output_arrays["maxresfield"]) = maxfield;
}

void OutputMatrices::setup_EH_and_gridlabels(SimulationParameters params, GridLabels input_grid_labels, PreferredInterpolationMethods pim) {
  // set the interpolation methods
  set_interpolation_methods(pim);

  // the dimensions of the fields
  E.il = H.il = (params.pml.Dxl) ? params.pml.Dxl + 2 : 0;
  E.iu = H.iu = (params.pml.Dxu) ? n_Yee_cells.I_tot() - params.pml.Dxu - 1 : n_Yee_cells.I_tot();
  E.jl = H.jl = (params.pml.Dyl) ? params.pml.Dyl + 2 : 0;
  E.ju = H.ju = (params.pml.Dyu) ? n_Yee_cells.J_tot() - params.pml.Dyu - 1 : n_Yee_cells.J_tot();
  E.kl = H.kl = (params.pml.Dzl) ? params.pml.Dzl + 2 : 0;
  E.ku = H.ku = (params.pml.Dzu) ? n_Yee_cells.K_tot() - params.pml.Dzu - 1 : n_Yee_cells.K_tot();
  // upper/lower limits of cell extraction
  E.I_tot = H.I_tot = E.iu - E.il + 1;
  E.J_tot = H.J_tot = E.ju - E.jl + 1;
  E.K_tot = H.K_tot = E.ku - E.kl + 1;

  // depending on the simulation, we might not need to assign this memory
  bool need_EH_memory = (params.run_mode == RunMode::complete && params.exphasorsvolume);
  // field matrices we will be assigning
  vector<string> field_matrices = {"Ex_out", "Ey_out", "Ez_out", "Hx_out", "Hy_out", "Hz_out"};
  // gridlabels matrices we will be assigning
  vector<string> grid_matrices = {"x_out", "y_out", "z_out"};
  if (!need_EH_memory) {
    // assign empty matrices
    output_arrays.assign_empty_matrix(field_matrices);
    output_arrays.assign_empty_matrix(grid_matrices);
  } else {
    // avoid memory leaks
    output_arrays.error_on_memory_assigned(field_matrices);
    output_arrays.error_on_memory_assigned(grid_matrices);

    // initialise to actual, proper arrays
    int dims[3] = {E.I_tot, E.J_tot, E.K_tot};
    spdlog::info("dims: ({0:d},{1:d},{2:d})", dims[0], dims[1], dims[2]);

    // create MATLAB data storage for the field-component outputs
    for(string matrix : field_matrices) {
      output_arrays[matrix] = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    }
    // create MATLAB data storage for the gridlabel outputs
    int label_dims_x[2] = {1, E.I_tot};
    int label_dims_y[2] = {1, E.J_tot};
    int label_dims_z[2] = {1, E.K_tot};
    output_arrays["x_out"] = mxCreateNumericArray(2, (const mwSize *) label_dims_x, mxDOUBLE_CLASS, mxREAL);
    output_arrays["y_out"] = mxCreateNumericArray(2, (const mwSize *) label_dims_y, mxDOUBLE_CLASS, mxREAL);
    output_arrays["z_out"] = mxCreateNumericArray(2, (const mwSize *) label_dims_z, mxDOUBLE_CLASS, mxREAL);

    // cast E, H, and gridlabel members to the MATLAB structures
    E.real.x = cast_matlab_3D_array(mxGetPr(output_arrays["Ex_out"]), E.I_tot, E.J_tot, E.K_tot);
    E.imag.x = cast_matlab_3D_array(mxGetPi(output_arrays["Ex_out"]), E.I_tot, E.J_tot, E.K_tot);

    E.real.y = cast_matlab_3D_array(mxGetPr(output_arrays["Ey_out"]), E.I_tot, E.J_tot, E.K_tot);
    E.imag.y = cast_matlab_3D_array(mxGetPi(output_arrays["Ey_out"]), E.I_tot, E.J_tot, E.K_tot);

    E.real.z = cast_matlab_3D_array(mxGetPr(output_arrays["Ez_out"]), E.I_tot, E.J_tot, E.K_tot);
    E.imag.z = cast_matlab_3D_array(mxGetPi(output_arrays["Ez_out"]), E.I_tot, E.J_tot, E.K_tot);

    H.real.x = cast_matlab_3D_array(mxGetPr(output_arrays["Hx_out"]), H.I_tot, H.J_tot, H.K_tot);
    H.imag.x = cast_matlab_3D_array(mxGetPi(output_arrays["Hx_out"]), H.I_tot, H.J_tot, H.K_tot);

    H.real.y = cast_matlab_3D_array(mxGetPr(output_arrays["Hy_out"]), H.I_tot, H.J_tot, H.K_tot);
    H.imag.y = cast_matlab_3D_array(mxGetPi(output_arrays["Hy_out"]), H.I_tot, H.J_tot, H.K_tot);

    H.real.z = cast_matlab_3D_array(mxGetPr(output_arrays["Hz_out"]), H.I_tot, H.J_tot, H.K_tot);
    H.imag.z = cast_matlab_3D_array(mxGetPi(output_arrays["Hz_out"]), H.I_tot, H.J_tot, H.K_tot);

    //now construct the grid labels
    output_grid_labels.x = mxGetPr(output_arrays["x_out"]);
    output_grid_labels.y = mxGetPr(output_arrays["y_out"]);
    output_grid_labels.z = mxGetPr(output_arrays["z_out"]);

    //initialise field arrays
    E.zero();
    H.zero();

    // initialise the gridlabels from the input grid labels (essentially shave off the PML labels)
    output_grid_labels.initialise_from(input_grid_labels, E.il, E.iu, E.jl, E.ju, E.kl, E.ku);
  }
}

void OutputMatrices::setup_Id(bool empty_allocation, int n_frequencies, int n_det_modes) {
  // create output by reserving memory
  if (empty_allocation) {
    output_arrays.assign_empty_matrix({"Id"});
  } else {
    // avoid memory leaks
    output_arrays.error_on_memory_assigned({"Id"});

    // create a structure array with two fields
    int dims[2] = {1, 1};
    const char *fieldnames[] = {"Idx", "Idy"};
    output_arrays["Id"] = mxCreateStructArray(2, (const mwSize *) dims, 2, fieldnames);

    // setup the IDVariables member so we can later write data to it
    ID.link_to_pointer(output_arrays["Id"], n_frequencies, n_det_modes);
  }
}

void OutputMatrices::setup_interpolation_outputs(SimulationParameters params) {
  // if we do not need to interpolate, we do not need to allocate memory to the interpolation-related outputs
  bool need_to_interpolate = (params.run_mode == RunMode::complete && params.exphasorsvolume);

  // output names that we will be allocating here
  vector<string> interpolated_gridlabels_names = {"x_i", "y_i", "z_i"};
  vector<string> interpolated_E_field_names = {"Ex_i", "Ey_i", "Ez_i"};
  vector<string> interpolated_H_field_names = {"Hx_i", "Hy_i", "Hz_i"};

  if (!need_to_interpolate) {
    // assign empty matrices
    output_arrays.assign_empty_matrix(interpolated_E_field_names);
    output_arrays.assign_empty_matrix(interpolated_H_field_names);
    output_arrays.assign_empty_matrix(interpolated_gridlabels_names);
  } else {
    // avoid memory leaks
    output_arrays.error_on_memory_assigned(interpolated_E_field_names);
    output_arrays.error_on_memory_assigned(interpolated_H_field_names);
    output_arrays.error_on_memory_assigned(interpolated_gridlabels_names);

    // interpolated field memory needs to be assigned. Since E, H field have the same {IJK}_tot, we can set the dimensions once, and do it here rather than separately for each field
    int i_upper = E.I_tot - 2, i_lower = 2;
    int j_upper = E.J_tot - 2, j_lower = 2;
    int k_upper = E.K_tot - 2, k_lower = 2;
    if (params.dimension != THREE) {
      k_upper = k_lower = 0;
    }
    // dimensions of the interpolated fields
    int outdims[3] = {i_upper - i_lower + 1, j_upper - j_lower + 1, 0};
    if (params.dimension == THREE) {
      outdims[2] = k_upper - k_lower + 1;
      outdims[1] = max(1, j_upper - j_lower + 1);
    }

    // Create interpolated E-field memory
    for (string matrix : interpolated_E_field_names) {
      output_arrays[matrix] =
              mxCreateNumericArray(3, (mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
    }
    // Create interpolated H-field memory
    for (string matrix : interpolated_H_field_names) {
      output_arrays[matrix] =
              mxCreateNumericArray(3, (mwSize *) outdims, mxDOUBLE_CLASS, mxCOMPLEX);
    }

    // interpolated gridlabels memory needs to be assigned
    // x-labels are always a fixed size
    int label_dims[2] = {1, E.I_tot - 3};
    output_arrays["x_i"] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);
    // y-labels might recieve nJ < 0, in which case a manual adjustment is made
    label_dims[1] = max(E.J_tot - 3, 1);
    output_arrays["y_i"] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);
    // z-labels depend on the dimensionality of the simulation
    if (params.dimension == THREE) {
      label_dims[1] = E.K_tot - 3;
    } else {
      label_dims[1] = 1;
    }
    output_arrays["z_i"] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);

    // now run the post-processing interpolation
    compute_interpolated_fields(params.dimension);
  }
}

void OutputMatrices::assign_surface_phasor_outputs(bool empty_allocation, mxArray *mx_surface_facets) {
  // output names that we will be assigning to
  vector<string> phasor_matrices = {"vertices", "camplitudes", "facets"};

  if (empty_allocation) {
    output_arrays.assign_empty_matrix(phasor_matrices, mxDOUBLE_CLASS, mxREAL);
  } else {
    // avoid memory leaks
    output_arrays.error_on_memory_assigned(phasor_matrices);

    // pull the memory blocks from the phasors
    output_arrays["vertices"] = surface_phasors.get_vertex_list();
    output_arrays["camplitudes"] =
            surface_phasors.get_mx_surface_amplitudes();
    output_arrays["facets"] = mx_surface_facets;
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
    int mpv_out = matPutVariable(output_file, matrix_to_write.c_str(), output_arrays[matrix_to_write]);
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
