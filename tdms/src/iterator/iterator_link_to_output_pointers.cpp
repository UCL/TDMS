#include "iterator_loop_variables.h"

// For mxGetPr, mxGetPi
#include "matrix.h"
// For cast_matlab_2D_array
#include "matlabio.h"

using namespace std;

void Iterator_LoopVariables::link_fields_and_labels_to_output_pointers(mxArray *output_pointers[]) {
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

void Iterator_LoopVariables::link_id_to_output_pointers(mxArray *output_pointers[], Iterator_IntermediateMATLABVariables &iMVars) {
  if (params.exdetintegral && params.run_mode == RunMode::complete) {
    int ndims = 2;
    int dims[2] = {1, 1};
    const char *fieldnames[] = {"Idx", "Idy"};
    // create output structure array
    output_pointers[26] = mxCreateStructArray(ndims, (const mwSize *) dims, 2, fieldnames);

    // now construct the fields for the structure array
    dims[0] = D_tilde.num_det_modes();
    dims[1] = f_ex_vec.size();

    iMVars.mx_Idx = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    iMVars.Idx_re = cast_matlab_2D_array(mxGetPr(iMVars.mx_Idx), dims[0], dims[1]);
    iMVars.Idx_im = cast_matlab_2D_array(mxGetPi(iMVars.mx_Idx), dims[0], dims[1]);

    iMVars.mx_Idy = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    iMVars.Idy_re = cast_matlab_2D_array(mxGetPr(iMVars.mx_Idy), dims[0], dims[1]);
    iMVars.Idy_im = cast_matlab_2D_array(mxGetPi(iMVars.mx_Idy), dims[0], dims[1]);

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
        iMVars.Idx_re[ifx][im] = 0.;
        iMVars.Idx_im[ifx][im] = 0.;
        iMVars.Idy_re[ifx][im] = 0.;
        iMVars.Idy_im[ifx][im] = 0.;
      }
    }

    mxSetField(output_pointers[26], 0, "Idx", iMVars.mx_Idx);
    mxSetField(output_pointers[26], 0, "Idy", iMVars.mx_Idy);
  } else {
    int ndims = 2;
    int dims[2] = {0, 0};
    output_pointers[26] =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  }
}
