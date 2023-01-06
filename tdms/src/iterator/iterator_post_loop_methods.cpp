#include "iterator.h"

#include "mesh_base.h"

void Iterator::normalise_field_volumes() {
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    E.normalise_volume();
    H.normalise_volume();
  }
}

void Iterator::normalise_surface_phasors() {
  if (params.run_mode == RunMode::complete && params.exphasorssurface) {
    spdlog::info("Surface phasors:");
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      surface_phasors.normalise_surface(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
    }
  }
}

void Iterator::normalise_vertex_phasors() {
  if (params.run_mode == RunMode::complete && vertex_phasors.there_are_vertices_to_extract_at()) {
    spdlog::info("Vertex phasors:");
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      vertex_phasors.normalise_vertices(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
    }
  }
}

void Iterator::normalise_Id_arrays() {
  if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete &&
      params.exdetintegral) {
    for (int im = 0; im < D_tilde.num_det_modes(); im++) {
      for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
        Idx[ifx][im] = Idx[ifx][im] / E_norm[ifx];
        Idy[ifx][im] = Idy[ifx][im] / E_norm[ifx];

        Idx_re[ifx][im] = real(Idx[ifx][im]);
        Idx_im[ifx][im] = imag(Idx[ifx][im]);

        Idy_re[ifx][im] = real(Idy[ifx][im]);
        Idy_im[ifx][im] = imag(Idy[ifx][im]);
      }
    }
  }
}

void Iterator::interpolate_field_values(mxArray *output_matrices[]) {
  if (params.dimension == THREE) {
    E.interpolate_over_range(&output_matrices[13], &output_matrices[14], &output_matrices[15], 2,
                             E.I_tot - 2, 2, E.J_tot - 2, 2, E.K_tot - 2, Dimension::THREE);
    H.interpolate_over_range(&output_matrices[16], &output_matrices[17], &output_matrices[18], 2,
                             H.I_tot - 2, 2, H.J_tot - 2, 2, H.K_tot - 2, Dimension::THREE);
  } else {
    // either TE or TM, but interpolate_over_range will handle that for us. Only difference is the k_upper/lower values we pass...
    E.interpolate_over_range(&output_matrices[13], &output_matrices[14], &output_matrices[15], 2,
                             E.I_tot - 2, 2, E.J_tot - 2, 0, 0, params.dimension);
    H.interpolate_over_range(&output_matrices[16], &output_matrices[17], &output_matrices[18], 2,
                             H.I_tot - 2, 2, H.J_tot - 2, 0, 0, params.dimension);
  }
}

void Iterator::write_interpolated_gridlabels(mxArray *output_matrices[]) {
  GridLabels interp_output_grid_labels;
  interp_output_grid_labels.x = mxGetPr((mxArray *) output_matrices[19]);
  interp_output_grid_labels.y = mxGetPr((mxArray *) output_matrices[20]);
  interp_output_grid_labels.z = mxGetPr((mxArray *) output_matrices[21]);

  if (params.dimension == THREE) {
    interp_output_grid_labels.initialise_from(output_grid_labels, 2, E.I_tot - 2, 2, E.J_tot - 2, 2,
                                              E.K_tot - 2);
  } else {
    interp_output_grid_labels.initialise_from(output_grid_labels, 2, E.I_tot - 2, 2, E.J_tot - 2, 0,
                                              0);
  }
}

void Iterator::regenerate_mesh_for_facets(mxArray *output_matrices[]) {
  mxArray *dummy_vertex_list;//< Temporary MATLAB storage

  // regenerate the mesh, since we threw away the facet list before iterating
  if (J_tot == 0)
    conciseCreateBoundary(cuboid[0], cuboid[1], cuboid[4], cuboid[5], &dummy_vertex_list,
                          &output_matrices[24]);
  else
    conciseTriangulateCuboidSkip(cuboid[0], cuboid[1], cuboid[2], cuboid[3], cuboid[4], cuboid[5],
                                 params.spacing_stride, &dummy_vertex_list, &output_matrices[24]);

  mxDestroyArray(dummy_vertex_list);//< Clear MATLAB memory

  //now create and populate the vertex list
  surface_phasors.create_vertex_list(input_grid_labels);

  //assign outputs (output_matrices[24] is assigned directly in the calls above)
  output_matrices[22] = surface_phasors.get_vertex_list();
  output_matrices[23] = surface_phasors.get_mx_surface_amplitudes();
}
