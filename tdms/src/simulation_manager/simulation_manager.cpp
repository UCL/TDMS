#include "simulation_manager.h"

#include <omp.h>
#include <spdlog/spdlog.h>

#include "mesh_base.h"

using namespace std;
using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;

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

double SimulationManager::linear_ramp(double t) {
  double period = 1. / (inputs.params.omega_an / (2 * DCPI));
  if (t > period * ramp_width) { return 1.; }
  else { return t / (period * ramp_width); }
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

SimulationManager::SimulationManager(InputMatrices in_matrices, SolverMethod _solver_method,
                                     PreferredInterpolationMethods _pim)
    : inputs(in_matrices, _solver_method, _pim) {
  solver_method = _solver_method;
  pim = _pim;

  // read number of Yee cells
  IJKDims IJK_tot = n_Yee_cells();

  // setup PSTD variables, and any dependencies there might be
  PSTD.set_using_dimensions(IJK_tot);
  if (solver_method == SolverMethod::PseudoSpectral) {
    int max_IJK = IJK_tot.max_IJK(), n_threads = omp_get_max_threads();
    eh_vec.allocate(n_threads, max_IJK + 1);

    inputs.E_s.initialise_fftw_plan(n_threads, eh_vec);
    inputs.H_s.initialise_fftw_plan(n_threads, eh_vec);
  }

  // these are needed later... but don't seem to EVER be used? They were previously plhs[6->9], but these outputs were never written. Also, they are assigned to, but never written out nor referrenced by any of the other variables in the main loop. I am confused... Also note that because we're using the Matrix class, we order indices [i][j][k] rather than [k][j][i] like in the rest of the codebase :(
  FDTD.allocate_memory(IJK_tot);

  // setup the output object
  prepare_output(in_matrices["fieldsample"], in_matrices["campssample"]);

  // initialise the {E,H}_norm variables to an array of zeros
  E_norm = vector<complex<double>>(inputs.f_ex_vec.size(), 0);
  H_norm = vector<complex<double>>(inputs.f_ex_vec.size(), 0);
}

void SimulationManager::post_loop_processing() {
  // normalise output fields
  if (inputs.params.run_mode == RunMode::complete && inputs.params.exphasorsvolume) {
    outputs.E.normalise_volume();
    outputs.H.normalise_volume();
  }

  // normalise the surface phasors
  if (inputs.params.run_mode == RunMode::complete && inputs.params.exphasorssurface) {
    spdlog::info("Surface phasors");
    for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
      outputs.surface_phasors.normalise_surface(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]),
                   imag(E_norm[ifx]));
    }
  }
  // normalise the vertex phasors
  if (inputs.params.run_mode == RunMode::complete &&
      outputs.vertex_phasors.there_are_vertices_to_extract_at()) {
    spdlog::info("Vertex phasors");
    for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
      outputs.vertex_phasors.normalise_vertices(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]),
                   imag(E_norm[ifx]));
    }
  }

  // write the ID output
  if (inputs.params.source_mode == SourceMode::pulsed &&
      inputs.params.run_mode == RunMode::complete && inputs.params.exdetintegral) {
    for (int im = 0; im < inputs.D_tilde.num_det_modes(); im++)
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        outputs.ID.x[ifx][im] = outputs.ID.x[ifx][im] / E_norm[ifx];
        outputs.ID.y[ifx][im] = outputs.ID.y[ifx][im] / E_norm[ifx];

        outputs.ID.x_real[ifx][im] = real(outputs.ID.x[ifx][im]);
        outputs.ID.x_imag[ifx][im] = imag(outputs.ID.x[ifx][im]);

        outputs.ID.y_real[ifx][im] = real(outputs.ID.y[ifx][im]);
        outputs.ID.y_imag[ifx][im] = imag(outputs.ID.y[ifx][im]);
      }
  }

  /* Find the maximum absolute value of residual field in the grid, and write too the outputs */
  double maxfield = std::max(inputs.E_s.largest_field_value(), inputs.H_s.largest_field_value());
  outputs.set_maxresfield(maxfield, false);

  // setup interpolated field outputs and labels (if necessary)
  outputs.setup_interpolation_outputs(inputs.params);

  /*Now export 3 matrices, a vertex list, a matrix of complex amplitudes at
    these vertices and a list of facets*/
  // if we need to extract the phasors, we will need to allocate memory in the output
  bool extracting_phasors =
          (inputs.params.exphasorssurface && inputs.params.run_mode == RunMode::complete);
  if (extracting_phasors) {
    //first regenerate the mesh since we threw away the facet list before iterating
    mxArray *dummy_vertex_list, *mx_surface_facets;
    if (n_Yee_cells().J_tot() == 0) {
      conciseCreateBoundary(inputs.cuboid[0], inputs.cuboid[1], inputs.cuboid[4], inputs.cuboid[5],
                            &dummy_vertex_list, &mx_surface_facets);
    } else {
      conciseTriangulateCuboidSkip(inputs.cuboid[0], inputs.cuboid[1], inputs.cuboid[2],
                                   inputs.cuboid[3], inputs.cuboid[4], inputs.cuboid[5],
                                   inputs.params.spacing_stride, &dummy_vertex_list,
                                   &mx_surface_facets);
    }
    // no longer need the dummy_vertex_list data
    mxDestroyArray(dummy_vertex_list);

    //now create and populate the vertex list
    outputs.surface_phasors.create_vertex_list(inputs.input_grid_labels);
    // assign to the output
    outputs.assign_surface_phasor_outputs(!extracting_phasors, mx_surface_facets);
    // now safe to cast mx_surface_vertices to nullptr if we need to, since surface_phasors will destroy the MATLAB memory
  }
}
