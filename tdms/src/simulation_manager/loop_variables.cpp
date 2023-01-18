#include "loop_variables.h"

#include <spdlog/spdlog.h>
#include <omp.h>

#include "globals.h"

using namespace tdms_phys_constants;
using namespace std;

LoopVariables::LoopVariables(ObjectsFromInfile &data, IJKDims E_field_dims, SolverMethod solver_method) {
  // deduce the number of non-pml cells in the z-direction, for efficiency
  K = data.IJK_tot.K_tot() - data.params.pml.Dxl - data.params.pml.Dxu;

  // deduce refractive index, and print to log
  refind = sqrt(1. / (data.freespace_Cbx[0] / data.params.dt * data.params.delta.dx) / EPSILON0);
  spdlog::info("refind (Refractive index) = {0:e}", refind);

  // setup temporary storgage for detector sensitivity evaluation
  if (data.params.exdetintegral) {
    int n0 = data.IJK_tot.I_tot() - data.params.pml.Dxl - data.params.pml.Dxu;
    int n1 = data.IJK_tot.J_tot() - data.params.pml.Dyl - data.params.pml.Dyu;
    Ex_t.initialise(n1, n0);
    Ey_t.initialise(n1, n0);
  }

  // We need to test for convergence under the following conditions. As such, we need to initialise the array (E_copy) that will ultimately be copies of the phasors at the previous iteration, to test convergence against
  if (data.params.run_mode == RunMode::complete && data.params.exphasorsvolume &&
      data.params.source_mode == SourceMode::steadystate) {
    int dummy_dims[3] = {E_field_dims.I_tot(), E_field_dims.J_tot(), E_field_dims.K_tot()};
    // allocate memory space for this array
    E_copy_MATLAB_data[0] =
            mxCreateNumericArray(3, (const mwSize *) dummy_dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
    E_copy_MATLAB_data[1] =
            mxCreateNumericArray(3, (const mwSize *) dummy_dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
    E_copy_MATLAB_data[2] =
            mxCreateNumericArray(3, (const mwSize *) dummy_dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez

    E_copy.real.x = cast_matlab_3D_array(mxGetPr(E_copy_MATLAB_data[0]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);
    E_copy.imag.x = cast_matlab_3D_array(mxGetPi(E_copy_MATLAB_data[0]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);

    E_copy.real.y = cast_matlab_3D_array(mxGetPr(E_copy_MATLAB_data[1]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);
    E_copy.imag.y = cast_matlab_3D_array(mxGetPi(E_copy_MATLAB_data[1]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);

    E_copy.real.z = cast_matlab_3D_array(mxGetPr(E_copy_MATLAB_data[2]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);
    E_copy.imag.z = cast_matlab_3D_array(mxGetPi(E_copy_MATLAB_data[2]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);

    E_copy.I_tot = E_field_dims.I_tot();
    E_copy.J_tot = E_field_dims.J_tot();
    E_copy.K_tot = E_field_dims.K_tot();

    E_copy.zero();
  }

  // Setup dispersive properties and the related arrays
  setup_dispersive_properties(data);

  // attempt to optimise the iteration loops, if we can
  optimise_loop_J_range(data);
}

void LoopVariables::setup_dispersive_properties(ObjectsFromInfile &data) {
  IJKDims IJK = data.IJK_tot;
  // determine whether or not we have a dispersive medium
  is_disp = is_dispersive(data.materials, IJK, data.gamma, data.params.dt);
  // work out if we have conductive background: background is conductive if at least one entry exceeds 1e-15
  is_conductive = !(data.rho_cond.all_elements_less_than(1e-15, IJK.I_tot() + 1, IJK.J_tot() + 1, IJK.K_tot() + 1));

  // prepare additional field variables for dispersive media
  E_nm1 = ElectricSplitField(IJK.I_tot(), IJK.J_tot(), IJK.K_tot());
  J_nm1 = CurrentDensitySplitField(IJK.I_tot(), IJK.J_tot(), IJK.K_tot());
  J_c = CurrentDensitySplitField(IJK.I_tot(), IJK.J_tot(), IJK.K_tot());
  // if we have a dispersive material we will need to write to the additional fields, so assign the memory to them and zero the entries
  if (is_disp || data.params.is_disp_ml) {
    E_nm1.allocate_and_zero();
    J_nm1.allocate_and_zero();
    J_s.allocate_and_zero();
  }
  // if we have a conductive material we will also need the conductivity/current-density of each cell
  if (is_conductive) { J_c.allocate_and_zero(); }
}

bool LoopVariables::is_dispersive(uint8_t ***materials, IJKDims IJK_tot,
                                  double *attenuation_constants, double dt, double non_zero_tol) {
  int max_mat = 0;
  // determine the number of entries in gamma, by examining the materials array
  for (int k = 0; k < (IJK_tot.K_tot() + 1); k++)
    for (int j = 0; j < (IJK_tot.J_tot() + 1); j++)
      for (int i = 0; i < (IJK_tot.I_tot() + 1); i++) {
        if (materials[k][j][i] > max_mat) max_mat = materials[k][j][i];
      }
  // now see if there are any non-zero attenuation constants
  for (int i = 0; i < max_mat; i++) {
    if (fabs(attenuation_constants[i] / dt) > non_zero_tol) { return true; }
  }
  return false;
}

LoopVariables::~LoopVariables() {
  for(int i = 0; i < 3; i++) {
    if (E_copy_MATLAB_data[i] != nullptr) { mxDestroyArray(E_copy_MATLAB_data[i]); }
  }
}

void LoopVariables::optimise_loop_J_range(ObjectsFromInfile &data, double non_zero_tol) {
  bool ksource_nz[4];//< "k source non-zero"
  for (int icomp = 0; icomp < 4; icomp++) { ksource_nz[icomp] = false; }

  if (data.IJK_tot.J_tot() == 0) {
    for (int icomp = 0; icomp < 4; icomp++)
      for (int ki = 0; ki < (data.IJK_tot.I_tot() + 1); ki++) {
        ksource_nz[icomp] = ksource_nz[icomp] ||
                            (fabs(data.Ksource.imag[0][ki - (data.I0.index)][icomp]) > non_zero_tol) ||
                            (fabs(data.Ksource.real[0][ki - (data.I0.index)][icomp]) > non_zero_tol);
      }
  }
  /* We now know the following information:
    Ey and Hx receive an input from the source condition only if ksource_nz[2] or ksource_nz[1]
    are non-zero.
    Ex and Hy receive an input from the source condition only if ksource_nz[3] or ksource_nz[0]
    are non-zero.

    Thus we can set the variables J_tot_p1_bound, J_tot_bound accordingly for the 3D, 2D-TE, and 2D-TM simulations.
  */
  J_tot_bound = data.IJK_tot.J_tot();
  J_tot_p1_bound = data.IJK_tot.J_tot() + 1;
  // If in a 2D simulation, adjust the values accordingly
  if (data.IJK_tot.J_tot() == 0) {
    // TE case
    if (ksource_nz[2] || ksource_nz[1] || data.params.eyi_present) {
      J_tot_bound = 1;
    } else {
      J_tot_bound = 0;
    }
    // TM case
    if (ksource_nz[3] || ksource_nz[0] || data.params.exi_present) {
      J_tot_p1_bound = 1;
    } else {
      J_tot_p1_bound = 0;
    }
  }
}
