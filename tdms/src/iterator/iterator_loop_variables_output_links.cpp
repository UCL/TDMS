#include "iterator_loop_variables.h"

#include <spdlog/spdlog.h>

// For mxGetPr, mxGetPi
#include "matrix.h"
// For cast_matlab_2D_array
#include "matlabio.h"

using namespace std;

void Iterator_LoopVariables::link_fields_and_labels(mxArray *output_pointers[]) {
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    // dimensions of the output arrays for the field and gridlabels
    int ndims = 3;
    mwSize dims[3] = {E.I_tot, E.J_tot, E.K_tot};
    mwSize label_dims[2] = {1, 1};
    spdlog::info("dims: ({0:d},{1:d},{2:d})", dims[0], dims[1], dims[2]);

    // create MATLAB data storage for the outputs
    output_pointers[0] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
    output_pointers[1] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
    output_pointers[2] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez

    output_pointers[3] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hx
    output_pointers[4] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hy
    output_pointers[5] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hz
    // Cast the E and H field C++ arrays to the MATLAB structures just created
    E.real.x = cast_matlab_3D_array(mxGetPr((mxArray *) output_pointers[0]), E.I_tot, E.J_tot,
                                    E.K_tot);
    E.imag.x = cast_matlab_3D_array(mxGetPi((mxArray *) output_pointers[0]), E.I_tot, E.J_tot,
                                    E.K_tot);

    E.real.y = cast_matlab_3D_array(mxGetPr((mxArray *) output_pointers[1]), E.I_tot, E.J_tot,
                                    E.K_tot);
    E.imag.y = cast_matlab_3D_array(mxGetPi((mxArray *) output_pointers[1]), E.I_tot, E.J_tot,
                                    E.K_tot);

    E.real.z = cast_matlab_3D_array(mxGetPr((mxArray *) output_pointers[2]), E.I_tot, E.J_tot,
                                    E.K_tot);
    E.imag.z = cast_matlab_3D_array(mxGetPi((mxArray *) output_pointers[2]), E.I_tot, E.J_tot,
                                    E.K_tot);

    H.real.x = cast_matlab_3D_array(mxGetPr((mxArray *) output_pointers[3]), H.I_tot, H.J_tot,
                                    H.K_tot);
    H.imag.x = cast_matlab_3D_array(mxGetPi((mxArray *) output_pointers[3]), H.I_tot, H.J_tot,
                                    H.K_tot);

    H.real.y = cast_matlab_3D_array(mxGetPr((mxArray *) output_pointers[4]), H.I_tot, H.J_tot,
                                    H.K_tot);
    H.imag.y = cast_matlab_3D_array(mxGetPi((mxArray *) output_pointers[4]), H.I_tot, H.J_tot,
                                    H.K_tot);

    H.real.z = cast_matlab_3D_array(mxGetPr((mxArray *) output_pointers[5]), H.I_tot, H.J_tot,
                                    H.K_tot);
    H.imag.z = cast_matlab_3D_array(mxGetPi((mxArray *) output_pointers[5]), H.I_tot, H.J_tot,
                                    H.K_tot);

    // The E_copy_data_placeholders need to be populated too, since these will be used to test steady state convergence. They'll ultimately end up as copies of the data in E at the timestep PRIOR to reaching convergence; we don't need to return them as an outputs, but do need them in the main iteration loop, so we free this memory in the destructor.
    if (params.source_mode == SourceMode::steadystate) {
      E_copy_data_placeholders[0] =
              mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
      E_copy_data_placeholders[1] =
              mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
      E_copy_data_placeholders[2] =
              mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez

      E_copy.real.x = cast_matlab_3D_array(mxGetPr((mxArray *) E_copy_data_placeholders[0]),
                                           dims[0], dims[1], dims[2]);
      E_copy.imag.x = cast_matlab_3D_array(mxGetPi((mxArray *) E_copy_data_placeholders[0]),
                                           dims[0], dims[1], dims[2]);

      E_copy.real.y = cast_matlab_3D_array(mxGetPr((mxArray *) E_copy_data_placeholders[1]),
                                           dims[0], dims[1], dims[2]);
      E_copy.imag.y = cast_matlab_3D_array(mxGetPi((mxArray *) E_copy_data_placeholders[1]),
                                           dims[0], dims[1], dims[2]);

      E_copy.real.z = cast_matlab_3D_array(mxGetPr((mxArray *) E_copy_data_placeholders[2]),
                                           dims[0], dims[1], dims[2]);
      E_copy.imag.z = cast_matlab_3D_array(mxGetPi((mxArray *) E_copy_data_placeholders[2]),
                                           dims[0], dims[1], dims[2]);

      E_copy.I_tot = E.I_tot;
      E_copy.J_tot = E.J_tot;
      E_copy.K_tot = E.K_tot;
    }

    // zero the newly-created arrays
    zero_field_arrays();

    // now construct the grid labels
    label_dims[1] = dims[0];
    output_pointers[10] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//x
    output_grid_labels.x = mxGetPr((mxArray *) output_pointers[10]);

    label_dims[1] = dims[1];
    output_pointers[11] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//y
    output_grid_labels.y = mxGetPr((mxArray *) output_pointers[11]);

    label_dims[1] = dims[2];
    output_pointers[12] =
            mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//y
    output_grid_labels.z = mxGetPr((mxArray *) output_pointers[12]);
  } else {
    //initialise to empty matrices
    int ndims = 2;
    int dims[2] = {0, 0};

    // The Ex, Ey, Ez, Hx, Hy, Hz components
    output_pointers[0] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[1] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[2] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[3] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[4] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[5] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);

    // {x,y,z}_grid_labels_out
    output_pointers[10] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[11] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    output_pointers[12] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  }
}

void Iterator_LoopVariables::link_id(mxArray *output_pointers[]) {
  if (params.exdetintegral && params.run_mode == RunMode::complete) {
    int ndims = 2;
    int dims[2] = {1, 1};
    const char *fieldnames[] = {"Idx", "Idy"};
    // create output structure array
    output_pointers[26] = mxCreateStructArray(ndims, (const mwSize *) dims, 2, fieldnames);

    // now construct the fields for the structure array
    dims[0] = D_tilde.num_det_modes();
    dims[1] = f_ex_vec.size();

    mx_Idx = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    Idx_re = cast_matlab_2D_array(mxGetPr(mx_Idx), dims[0], dims[1]);
    Idx_im = cast_matlab_2D_array(mxGetPi(mx_Idx), dims[0], dims[1]);

    mx_Idy = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    Idy_re = cast_matlab_2D_array(mxGetPr(mx_Idy), dims[0], dims[1]);
    Idy_im = cast_matlab_2D_array(mxGetPi(mx_Idy), dims[0], dims[1]);

    Idx = (complex<double> **) malloc(sizeof(complex<double> *) * f_ex_vec.size());
    Idy = (complex<double> **) malloc(sizeof(complex<double> *) * f_ex_vec.size());

    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      Idx[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
      Idy[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
      for (int im = 0; im < dims[0]; im++) {
        Idx[ifx][im] = 0.;
        Idy[ifx][im] = 0.;
      }
    }

    for (int im = 0; im < dims[0]; im++) {
      for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
        Idx_re[ifx][im] = 0.;
        Idx_im[ifx][im] = 0.;
        Idy_re[ifx][im] = 0.;
        Idy_im[ifx][im] = 0.;
      }
    }

    mxSetField(output_pointers[26], 0, "Idx", mx_Idx);
    mxSetField(output_pointers[26], 0, "Idy", mx_Idy);
  } else {
    int ndims = 2;
    int dims[2] = {0, 0};
    output_pointers[26] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  }
}

void Iterator_LoopVariables::link_fdtd_phasor_arrays(mxArray *output_pointers[]) {
  const int ndims = 2;
  int dims_ex_hy[2] = {I_tot, J_tot + 1};
  int dims_ey_hx[2] = {I_tot + 1, J_tot};

  // x electric field source phasor - boot strapping
  output_pointers[6] =
          mxCreateNumericArray(ndims, (const mwSize *) dims_ex_hy, mxDOUBLE_CLASS, mxCOMPLEX);
  iwave_lEx_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) output_pointers[6]), dims_ex_hy[0],
                                       dims_ex_hy[1]);
  iwave_lEx_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) output_pointers[6]), dims_ex_hy[0],
                                       dims_ex_hy[1]);
  zero_cast_array(iwave_lEx_Rbs, dims_ex_hy[0], dims_ex_hy[1]);
  zero_cast_array(iwave_lEx_Ibs, dims_ex_hy[0], dims_ex_hy[1]);

  // y electric field source phasor - boot strapping
  output_pointers[7] =
          mxCreateNumericArray(ndims, (const mwSize *) dims_ey_hx, mxDOUBLE_CLASS, mxCOMPLEX);
  iwave_lEy_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) output_pointers[7]), dims_ey_hx[0],
                                       dims_ey_hx[1]);
  iwave_lEy_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) output_pointers[7]), dims_ey_hx[0],
                                       dims_ey_hx[1]);
  zero_cast_array(iwave_lEy_Rbs, dims_ey_hx[0], dims_ey_hx[1]);
  zero_cast_array(iwave_lEy_Ibs, dims_ey_hx[0], dims_ey_hx[1]);

  // x magnetic field source phasor - boot strapping
  output_pointers[8] =
          mxCreateNumericArray(ndims, (const mwSize *) dims_ey_hx, mxDOUBLE_CLASS, mxCOMPLEX);

  iwave_lHx_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) output_pointers[8]), dims_ey_hx[0],
                                       dims_ey_hx[1]);
  iwave_lHx_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) output_pointers[8]), dims_ey_hx[0],
                                       dims_ey_hx[1]);
  zero_cast_array(iwave_lHx_Rbs, dims_ey_hx[0], dims_ey_hx[1]);
  zero_cast_array(iwave_lHx_Ibs, dims_ey_hx[0], dims_ey_hx[1]);

  // y magnetic field source phasor - boot strapping
  output_pointers[9] =
          mxCreateNumericArray(ndims, (const mwSize *) dims_ex_hy, mxDOUBLE_CLASS, mxCOMPLEX);
  iwave_lHy_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) output_pointers[9]), dims_ex_hy[0],
                                       dims_ex_hy[1]);
  iwave_lHy_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) output_pointers[9]), dims_ex_hy[0],
                                       dims_ex_hy[1]);
  zero_cast_array(iwave_lHy_Rbs, dims_ex_hy[0], dims_ex_hy[1]);
  zero_cast_array(iwave_lHy_Ibs, dims_ex_hy[0], dims_ex_hy[1]);
}

void Iterator_LoopVariables::link_fieldsample(mxArray *output_pointers[]) {
  output_pointers[27] = fieldsample.mx;
}

void Iterator_LoopVariables::link_vertex_phasors(mxArray *output_pointers[]) {
  vertex_phasors.setup_camplitude_arrays(f_ex_vec.size());
  output_pointers[28] = vertex_phasors.get_mx_camplitudes();
}
