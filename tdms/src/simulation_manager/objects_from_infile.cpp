#include "simulation_manager/objects_from_infile.h"

#include <spdlog/spdlog.h>
#include <stdexcept>

// For ptr_to_vector_in, ptr_to_vector_or_empty_in, int_cast_from_double_in
#include "matlabio.h"
// for init_grid_arrays
#include "array_init.h"

using tdms_math_constants::DCPI;

IndependentObjectsFromInfile::IndependentObjectsFromInfile(InputMatrices matrices_from_input_file,
                                                           SolverMethod _solver_method,
                                                           PreferredInterpolationMethods _pim)
    :// initialisation list - members whose classes have no default constructors
      Cmaterial(matrices_from_input_file["Cmaterial"]),// get Cmaterial
      Dmaterial(matrices_from_input_file["Dmaterial"]),// get Dmaterial
      C(matrices_from_input_file["C"]),                // get C
      D(matrices_from_input_file["D"]),                // get D
      I0(matrices_from_input_file["interface"], "I0"), // get the interface(s)
      I1(matrices_from_input_file["interface"], "I1"),
      J0(matrices_from_input_file["interface"], "J0"),
      J1(matrices_from_input_file["interface"], "J1"),
      K0(matrices_from_input_file["interface"], "K0"),
      K1(matrices_from_input_file["interface"], "K1"),
      matched_layer(matrices_from_input_file["dispersive_aux"]),// get dispersive_aux
      Ei(matrices_from_input_file["tdfield"])                   // get tdfield
{
  // set solver method
  set_solver_method(_solver_method);
  // set interpolation methods
  set_interpolation_method(_pim);

  // unpack the parameters for this simulation
  params.unpack_from_input_matrices(matrices_from_input_file);

  // get the fdtd grid
  init_grid_arrays(matrices_from_input_file["fdtdgrid"], E_s, H_s, materials);
  // set the {IJK}_tot variables using the split-field information we just unpacked
  IJK_tot = IJKDims(E_s.I_tot, E_s.J_tot, E_s.K_tot);

  // Get freespace - Cby Cbz Dbx Dby Dbz are unused
  freespace_Cbx =
          mxGetPr(ptr_to_vector_in(matrices_from_input_file["freespace"], "Cbx", "freespace"));

  // Get disp_params
  alpha = mxGetPr(ptr_to_vector_or_empty_in(matrices_from_input_file["disp_params"], "alpha",
                                            "disp_params"));
  beta = mxGetPr(ptr_to_vector_or_empty_in(matrices_from_input_file["disp_params"], "beta",
                                           "disp_params"));
  gamma = mxGetPr(ptr_to_vector_or_empty_in(matrices_from_input_file["disp_params"], "gamma",
                                            "disp_params"));

  // Get grid_labels (not needed in main loop)
  input_grid_labels = GridLabels(matrices_from_input_file["grid_labels"]);

  // Get phasorsurface
  cuboid = Cuboid();
  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    cuboid.initialise(matrices_from_input_file["phasorsurface"], IJK_tot.J_tot());
  }

  // Get conductive_aux, and setup with pointers
  // mxGetPr pointers will be cleaned up by XYZVectors destructor
  rho_cond = XYZVectors();
  rho_cond.x = mxGetPr(
          ptr_to_vector_in(matrices_from_input_file["conductive_aux"], "rho_x", "conductive_aux"));
  rho_cond.y = mxGetPr(
          ptr_to_vector_in(matrices_from_input_file["conductive_aux"], "rho_y", "conductive_aux"));
  rho_cond.z = mxGetPr(
          ptr_to_vector_in(matrices_from_input_file["conductive_aux"], "rho_z", "conductive_aux"));

  // prepare variables dependent on frequency extraction vector
  f_vec = FrequencyVectors();
  pupil = Pupil();
  D_tilde = DTilde();
  // if exdetintegral is flagged, setup pupil, D_tilde, and f_vec accordingly
  if (params.exdetintegral) {
    f_vec.initialise(matrices_from_input_file["f_vec"]);
    pupil.initialise(matrices_from_input_file["Pupil"], f_vec.x.size(), f_vec.y.size());
    D_tilde.initialise(matrices_from_input_file["D_tilde"], f_vec.x.size(), f_vec.y.size());

    if (!mxIsEmpty(matrices_from_input_file["k_det_obs_global"])) {
      params.k_det_obs =
              int_cast_from_double_in(matrices_from_input_file["k_det_obs_global"], "k_det_obs") -
              1;
    }
    params.z_obs = input_grid_labels.z[params.k_det_obs];
  }

  // Get tdfdir - check what's going on here in a second
  ex_td_field_exporter = TDFieldExporter2D();
  if (mxIsChar(matrices_from_input_file["tdfdir"])) {
    ex_td_field_exporter.folder_name =
            string_in(matrices_from_input_file["tdfdir"], "tdfdir").c_str();

    int Ni_tdf = 0, Nk_tdf = 0;
    for (int k = 0; k < IJK_tot.K_tot(); k++)
      if ((k % skip_tdf) == 0) Nk_tdf++;

    for (int i = 0; i < IJK_tot.I_tot(); i++)
      if ((i % skip_tdf) == 0) Ni_tdf++;
    spdlog::info("Ni_tdf = {0:d}, Nk_tdf = {1:d}", Ni_tdf, Nk_tdf);

    if (!are_equal(ex_td_field_exporter.folder_name, "")) {
      params.has_tdfdir = true;
      ex_td_field_exporter.allocate(Ni_tdf, Nk_tdf);
    }
  }

  // work out if we have a dispersive background
  if (params.is_disp_ml) { params.is_disp_ml = matched_layer.is_dispersive(IJK_tot.K_tot()); }

  // Set dt so that an integer number of time periods fits within a sinusoidal period
  double Nsteps_tmp = 0.0;
  double dt_old;
  if (params.source_mode == SourceMode::steadystate) {
    dt_old = params.dt;
    Nsteps_tmp = ceil(2. * DCPI / params.omega_an / params.dt * 3);
    params.dt = 2. * DCPI / params.omega_an * 3 / Nsteps_tmp;
    if (params.run_mode == RunMode::complete) {
      spdlog::info("Changed dt to {0:.10e} (was {1:.10e})", params.dt, dt_old);
    }
  }
  Nsteps = (int) lround(Nsteps_tmp);
  // Nt should be an integer number of Nsteps in the case of steady-state operation
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete) {
    if (params.Nt / Nsteps * Nsteps != params.Nt) {
      int old_Nt = params.Nt;//< For logging purposes, holds the Nt value that had to be changed
      params.Nt = params.Nt / Nsteps * Nsteps;
      spdlog::info("Changing the value of Nt to {0:d} (was {1:d})", params.Nt, old_Nt);
    }
    spdlog::info("Nsteps: {0:d}", Nsteps);
  }
}

IndependentObjectsFromInfile::~IndependentObjectsFromInfile() {
  if (params.dimension == THREE) {
    free_cast_matlab_3D_array(materials, E_s.K_tot + 1);
  } else {
    free_cast_matlab_3D_array(materials, 0);
  }
}

ObjectsFromInfile::ObjectsFromInfile(InputMatrices matrices_from_input_file,
                                     SolverMethod _solver_method,
                                     PreferredInterpolationMethods _pim)
    :// build the independent objects first
      IndependentObjectsFromInfile(matrices_from_input_file, _solver_method, _pim),
      // Source has no default constructor, and we need information from the Iterator_IndependentObjectsFromInfile first
      Isource(matrices_from_input_file["Isource"], J1.index - J0.index + 1, K1.index - K0.index + 1,
              "Isource"),
      Jsource(matrices_from_input_file["Jsource"], I1.index - I0.index + 1, K1.index - K0.index + 1,
              "Jsource"),
      Ksource(matrices_from_input_file["Ksource"], I1.index - I0.index + 1, J1.index - J0.index + 1,
              "Ksource"),
      // Get structure, we need I_tot from Iterator_IndependentObjectsFromInfile
      structure(matrices_from_input_file["structure"], IJK_tot.I_tot()),
      // Get f_ex_vec, the vector of frequencies to extract the field at.
      // Need params.omega from Iterator_IndependentObjectsFromInfile
      f_ex_vec(matrices_from_input_file["f_ex_vec"], params.omega_an) {
  // Update params according to structure's values
  params.is_structure = structure.has_elements();
}

ObjectsFromInfile::~ObjectsFromInfile() {
  if (I0.apply || I1.apply) {
    free_cast_matlab_3D_array(Isource.imag, (K1.index - K0.index + 1));
    free_cast_matlab_3D_array(Isource.real, (K1.index - K0.index + 1));
  }
  if (J0.apply || J1.apply) {
    free_cast_matlab_3D_array(Jsource.imag, (K1.index - K0.index + 1));
    free_cast_matlab_3D_array(Jsource.real, (K1.index - K0.index + 1));
  }
  if (K0.apply || K1.apply) {
    free_cast_matlab_3D_array(Ksource.imag, (J1.index - J0.index + 1));
    free_cast_matlab_3D_array(Ksource.real, (J1.index - J0.index + 1));
  }
}
