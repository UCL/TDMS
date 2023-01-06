#include "iterator_objects_from_infile.h"

#include <stdexcept>

// For ptr_to_vector_in, ptr_to_vector_or_empty_in, int_cast_from_double_in
#include "matlabio.h"
// for init_grid_arrays
#include "array_init.h"

Iterator_IndependentObjectsFromInfile::Iterator_IndependentObjectsFromInfile(
        InputMatrices matrices_from_input_file, SolverMethod _solver_method)
    :// initalizer list - members whose classes have no default constructors
      Cmaterial(matrices_from_input_file["Cmaterial"]),// get Cmaterial
      Dmaterial(matrices_from_input_file["Dmaterial"]),// get Dmaterial
      C(matrices_from_input_file["C"]),                // get C
      D(matrices_from_input_file["D"]),                // get D
      I0(matrices_from_input_file["interface"], "I0"), // get the interface
      I1(matrices_from_input_file["interface"], "I1"),
      J0(matrices_from_input_file["interface"], "J0"),
      J1(matrices_from_input_file["interface"], "J1"),
      K0(matrices_from_input_file["interface"], "K0"),
      K1(matrices_from_input_file["interface"], "K1"),
      matched_layer(matrices_from_input_file["dispersive_aux"]),// get dispersive_aux
      Ei(matrices_from_input_file["tdfield"]),                  // get tdfield
      fieldsample(matrices_from_input_file["fieldsample"]),     // get fieldsample
      vertex_phasors(matrices_from_input_file["campssample"])   // get vertex_phasors
{
  // set skip_tdf based on the solver method
  solver_method = _solver_method;
  if (solver_method == SolverMethod::FiniteDifference) {
    skip_tdf = 6;
  } else if (solver_method == SolverMethod::PseudoSpectral) {
    skip_tdf = 1;
  } else {
    throw std::runtime_error("Solver method not recognised!");
  }

  // unpack the parameters for this simulation
  params.unpack_from_input_matrices(matrices_from_input_file);

  // get the fdtd grid
  init_grid_arrays(matrices_from_input_file["fdtdgrid"], E_s, H_s, materials);
  // set the {IJK}_tot variables using the split-field information we just unpacked
  set_IJK_tot_from_E_s();

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
    cuboid.initialise(matrices_from_input_file["phasorsurface"], J_tot);
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
    for (int k = 0; k < K_tot; k++)
      if ((k % skip_tdf) == 0) Nk_tdf++;

    for (int i = 0; i < I_tot; i++)
      if ((i % skip_tdf) == 0) Ni_tdf++;
    spdlog::info("Ni_tdf = {0:d}, Nk_tdf = {1:d}", Ni_tdf, Nk_tdf);

    if (!are_equal(ex_td_field_exporter.folder_name, "")) {
      params.has_tdfdir = true;
      ex_td_field_exporter.allocate(Ni_tdf, Nk_tdf);
    }
  }
}

Iterator_IndependentObjectsFromInfile::~Iterator_IndependentObjectsFromInfile() {
  spdlog::info("Destroying Iterator_IndependentObjectsFromInfile");

  if (params.dimension == THREE) {
    free_cast_matlab_3D_array(materials, E_s.K_tot + 1);
  } else {
    free_cast_matlab_3D_array(materials, 0);
  }
}

Iterator_ObjectsFromInfile::Iterator_ObjectsFromInfile(InputMatrices matrices_from_input_file,
                                                       SolverMethod _solver_method)
    : // build the independent objects first
      Iterator_IndependentObjectsFromInfile(matrices_from_input_file, _solver_method),
      // Source has no default constructor, and we need information from the Iterator_IndependentObjectsFromInfile first
      Isource(matrices_from_input_file["Isource"], J1.index - J0.index + 1, K1.index - K0.index + 1,
              "Isource"),
      Jsource(matrices_from_input_file["Jsource"], I1.index - I0.index + 1, K1.index - K0.index + 1,
              "Jsource"),
      Ksource(matrices_from_input_file["Ksource"], I1.index - I0.index + 1, J1.index - J0.index + 1,
              "Ksource"),
      // Get structure, we need I_tot from Iterator_IndependentObjectsFromInfile
      structure(matrices_from_input_file["structure"], I_tot),
      // Get f_ex_vec, the vector of frequencies to extract the field at.
      // Need params.omega from Iterator_IndependentObjectsFromInfile
      f_ex_vec(matrices_from_input_file["f_ex_vec"], params.omega_an)
{
  // Update params according to structure's values
  params.is_structure = structure.has_elements();
}

Iterator_ObjectsFromInfile::~Iterator_ObjectsFromInfile() {
  spdlog::info("Destroying Iterator_ObjectsFromInfile");
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
