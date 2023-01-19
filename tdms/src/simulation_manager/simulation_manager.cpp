#include "simulation_manager.h"

#include <omp.h>

using namespace std;
using namespace tdms_math_constants;

SimulationManager::SimulationManager(InputMatrices in_matrices, SolverMethod _solver_method,
                                     PreferredInterpolationMethods _pim)
    : inputs(in_matrices, _solver_method, _pim) {
  solver_method = _solver_method;
  pim = _pim;

  // read number of Yee cells
  IJKDims IJK_tot = n_Yee_cells();

  // setup PSTD variables, and any dependencies there might be
  if (solver_method == SolverMethod::PseudoSpectral) {
    int max_IJK = IJK_tot.max_IJK(), n_threads = omp_get_max_threads();

    PSTD.set_using_dimensions(IJK_tot);

    eh_vec.allocate(n_threads, max_IJK + 1);
    inputs.E_s.initialise_fftw_plan(n_threads, eh_vec);
    inputs.H_s.initialise_fftw_plan(n_threads, eh_vec);
  }

  // initialise the {E,H}_norm variables to an array of zeros
  E_norm = vector<complex<double>>(inputs.f_ex_vec.size(), 0);
  H_norm = vector<complex<double>>(inputs.f_ex_vec.size(), 0);

  // these are needed later... but don't seem to EVER be used? They were previously plhs[6->9], but these outputs were never written. Also, they are assigned to, but never written out nor referrenced by any of the other variables in the main loop. I am confused... Also note that because we're using the Matrix class, we order indices [i][j][k] rather than [k][j][i] like in the rest of the codebase :(
  FDTD = FDTDBootstrapper(IJK_tot);

  // setup the output object
  prepare_output(in_matrices["fieldsample"], in_matrices["campssample"]);
}

void SimulationManager::extract_phasor_norms(int frequency_index, int tind, int Nt) {
  double omega = inputs.f_ex_vec[frequency_index] * 2 * DCPI;
  E_norm[frequency_index] +=
          outputs.E.ft *
          exp(fmod(omega * ((double) (tind + 1)) * inputs.params.dt, 2 * DCPI) * IMAGINARY_UNIT) *
          1. / ((double) Nt);
  H_norm[frequency_index] +=
          outputs.H.ft *
          exp(fmod(omega * ((double) tind + 0.5) * inputs.params.dt, 2 * DCPI) * IMAGINARY_UNIT) *
          1. / ((double) Nt);
}

void SimulationManager::prepare_output(const mxArray *fieldsample, const mxArray *campssample) {
  IJKDims IJK_tot = n_Yee_cells();

  outputs.set_n_Yee_cells(IJK_tot);
  /*set up surface mesh if required*/
  outputs.setup_surface_mesh(inputs.cuboid, inputs.params, inputs.f_ex_vec.size());
  /*Now set up the phasor array, we will have 3 complex output arrays for Ex, Ey and Ez.
    Phasors are extracted over the range Dxl + 3 - 1 to I_tot - Dxu - 1 to avoid pml cells
    see page III.80 for explanation of the following. This has been extended so that interpolation
    is done at the end of the FDTD run and also to handle the case of when there is no PML in place
    more appropriatley*/
  outputs.setup_EH_and_gridlabels(inputs.params, inputs.input_grid_labels, pim);
  // Setup the ID output
  bool need_Id_memory =
          (inputs.params.exdetintegral && inputs.params.run_mode == RunMode::complete);
  outputs.setup_Id(!need_Id_memory, inputs.f_ex_vec.size(), inputs.D_tilde.num_det_modes());
  // Link to the fieldsample and vertex_phasor inputs
  outputs.setup_fieldsample(fieldsample);
  outputs.setup_vertex_phasors(campssample, inputs.f_ex_vec.size());
}
