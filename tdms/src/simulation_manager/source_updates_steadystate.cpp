#include "simulation_manager/simulation_manager.h"

#include "cell_coordinate.h"

using namespace std;
using namespace tdms_phys_constants;
using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;

void SimulationManager::update_source_terms_steadystate(
        double time_H, bool is_conductive, CurrentDensitySplitField &J_c,
        CurrentDensitySplitField &J_s) {
  update_Isource_terms_steadystate(time_H, is_conductive, J_c, J_s);
  update_Jsource_terms_steadystate(time_H, is_conductive, J_c, J_s);
  update_Ksource_terms_steadystate(time_H, is_conductive, J_c, J_s);
}

void SimulationManager::update_Isource_terms_steadystate(
        double time_H, bool is_conductive, CurrentDensitySplitField &J_c,
        CurrentDensitySplitField &J_s) {
  // Only run update equations is source data was provided
  if (inputs.Isource.is_empty()) { return; }

  int array_ind;

  // Update across I0, provided a source term is here
  if (inputs.I0.apply) {
    for (int k = (inputs.K0.index); k <= (inputs.K1.index); k++) {
      for (int j = (inputs.J0.index); j <= (inputs.J1.index); j++) {
        if (!inputs.params.is_multilayer) {
          array_ind = inputs.I0.index;
        } else {
          array_ind = (n_Yee_cells().i + 1) * k + inputs.I0.index;
        }
        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          update_source_steadystate(time_H, AxialDirection::X, true, true,
                                    is_conductive, j, k, array_ind, J_c, J_s);
        }
        if (j < (inputs.J1.index)) {
          update_source_steadystate(time_H, AxialDirection::X, false, true,
                                    is_conductive, j, k, array_ind, J_c, J_s);
        }
      }
    }
  }

  // Update across I1, provided a source term is here
  if (inputs.I1.apply) {
    for (int k = (inputs.K0.index); k <= (inputs.K1.index); k++) {
      for (int j = (inputs.J0.index); j <= (inputs.J1.index); j++) {
        if (!inputs.params.is_multilayer) {
          array_ind = inputs.I1.index;
        } else {
          array_ind = (n_Yee_cells().i + 1) * k + inputs.I1.index;
        }
        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          update_source_steadystate(time_H, AxialDirection::X, true, false,
                                    is_conductive, j, k, array_ind, J_c, J_s);
        }
        if (j < (inputs.J1.index)) {
          update_source_steadystate(time_H, AxialDirection::X, false, false,
                                    is_conductive, j, k, array_ind, J_c, J_s);
        }
      }
    }
  }
}

void SimulationManager::update_Jsource_terms_steadystate(
        double time_H, bool is_conductive, CurrentDensitySplitField &J_c,
        CurrentDensitySplitField &J_s) {
  // Only run update equations is source data was provided
  if (inputs.Jsource.is_empty()) { return; }

  int array_ind;

  // Update across J0, provided a source term is there
  if (inputs.J0.apply) {
    for (int k = inputs.K0.index; k <= inputs.K1.index; k++) {
      for (int i = inputs.I0.index; i < inputs.I1.index; i++) {
        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          // shouldn't this be one-level up to be consistent with the other
          // source update steps?!?!?!?!
          if (!inputs.params.is_multilayer) {
            array_ind = inputs.J0.index;
          } else {
            array_ind = (n_Yee_cells().j + 1) * k + inputs.J0.index;
          }

          update_source_steadystate(time_H, AxialDirection::Y, true, true,
                                    is_conductive, i, k, array_ind, J_c, J_s);
        }
        if (i < inputs.I1.index) {
          update_source_steadystate(time_H, AxialDirection::Y, false, true,
                                    is_conductive, i, k, array_ind, J_c, J_s);
        }
      }
    }
  }

  // Update across J1, provided a source term is there
  if (inputs.J1.apply) {
    for (int k = inputs.K0.index; k <= inputs.K1.index; k++) {
      for (int i = inputs.I0.index; i <= inputs.I1.index; i++) {
        if (!inputs.params.is_multilayer) {
          array_ind = inputs.J1.index;
        } else {
          array_ind = (n_Yee_cells().j + 1) * k + inputs.J1.index;
        }
        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          update_source_steadystate(time_H, AxialDirection::Y, true, false,
                                    is_conductive, i, k, array_ind, J_c, J_s);
        }
        if (i < (inputs.I1.index)) {
          update_source_steadystate(time_H, AxialDirection::Y, false, false,
                                    is_conductive, i, k, array_ind, J_c, J_s);
        }
      }
    }
  }
}

void SimulationManager::update_Ksource_terms_steadystate(
        double time_H, bool is_conductive, CurrentDensitySplitField &J_c,
        CurrentDensitySplitField &J_s) {
  // Only run update equations is source data was provided
  if (inputs.Ksource.is_empty()) { return; }

  // Update across K0, provided a source term is there
  if (inputs.K0.apply) {
    for (int j = inputs.J0.index; j <= inputs.J1.index; j++) {
      for (int i = inputs.I0.index; i <= inputs.I1.index; i++) {
        if (j < (inputs.J1.index)) {
          update_source_steadystate(time_H, AxialDirection::Z, true, true,
                                    is_conductive, i, j, inputs.K0.index, J_c,
                                    J_s);
        }
        if (i < (inputs.I1.index)) {
          update_source_steadystate(time_H, AxialDirection::Z, false, true,
                                    is_conductive, i, j, inputs.K0.index, J_c,
                                    J_s);
        }
      }
    }
  }

  // Update across K1, provided a source term is there
  if (inputs.K1.apply) {
    for (int j = inputs.J0.index; j <= inputs.J1.index; j++) {
      for (int i = inputs.I0.index; i <= inputs.I1.index; i++) {
        if (j < (inputs.J1.index)) {
          update_source_steadystate(time_H, AxialDirection::Z, true, false,
                                    is_conductive, i, j, inputs.K1.index, J_c,
                                    J_s);
        }
        if (i < (inputs.I1.index)) {
          update_source_steadystate(time_H, AxialDirection::Z, false, false,
                                    is_conductive, i, j, inputs.K1.index, J_c,
                                    J_s);
        }
      }
    }
  }
}

void SimulationManager::update_source_steadystate(
        double time_H, AxialDirection parallel, bool C_axis, bool zero_plane,
        bool is_conductive, int cell_b, int cell_c, int array_ind,
        CurrentDensitySplitField &J_c, CurrentDensitySplitField &J_s) {
  // Common phase and amplitude terms to apply in all computations
  complex<double> common_phase = exp(
          -IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_H, 2. * DCPI));
  double common_amplitude = linear_ramp(time_H);

  double c_constant;    //< C.b.{x,y,z} constant that appears in the update
                        // equations
  double conductive_aux;//< Conductivity, only present in a conductive material.
                        //= inputs.rho_cond.{x,y,z}

  /*! Dispersive parameter. = inputs.matched_layer.kappa.{x,y,z} */
  double kappa;
  /*! NOTE: if dealing with Ksource, then the index used to fetch the gamma
   * parameter is k = inputs.K1.index, rather than cell_c. This inconsistency &
   * potential error is flagged in #221. Have a feeling that, in the K-case,
   * this should be inputs.K{0,1}.index depending on the plane we're in */
  int gamma_ind = (parallel == AxialDirection::Z) ? inputs.K1.index : cell_c;
  /*! Dispersive parameter. = inputs.matched_layer.gamma[k]. */
  double gamma = inputs.matched_layer.gamma[gamma_ind];

  /*! Value of the source term in Yee cell (cell_a, cell_b, cell_c) */
  complex<double> source_value;
  /*! The split-field cell index to update */
  CellCoordinate cell_to_update;
  /*! Field updates are signed based on the component, plane, and axis. [0] =
   * E_s sign, [1] = J_c sign, [2] = J_s sign */
  double update_sign[3] = {1., 1., 1.};

  SplitFieldComponent
          *E_s_component,//< The E_s component that should be updated
          *J_c_component,//< The J_c component that should be updated
          *J_s_component;//< The J_s component that should be updated

  // If we are currently updating along the C-axis, the split_field_ID is one
  // less than if we are updating along the B-axis
  int split_field_ID = (C_axis) ? 2 : 3;
  // If we are updating in the {IJK}1, rather than the {IJK}0 plane, the
  // split_field_ID is 4 greater than usual (uses indices 6/7 as opposed to 2/3)
  if (!zero_plane) { split_field_ID += 4; }
  /*! This is the index that is to be passed to the (relevent) Source instance,
   * to recover the source value at the relevant Yee cell.

   cell_b and cell_c need to be offset by the corresponding {IJK}0.index before
   we cna retrieve the source value. */
  SourceIndex s_index = {split_field_ID, cell_b, cell_c};

  // setup axis-dependent parameters
  switch (parallel) {
    case AxialDirection::X:
      // Dealing with Isource. j = cell_b, k = cell_c
      c_constant = inputs.C.b.x[array_ind];

      s_index.j -= inputs.J0.index;
      s_index.k -= inputs.K0.index;
      source_value = inputs.Isource[s_index];

      cell_to_update = {inputs.I0.index, cell_b, cell_c};
      if (!zero_plane) { cell_to_update.i = inputs.I1.index; }

      conductive_aux = inputs.rho_cond.x[array_ind];
      kappa = inputs.matched_layer.kappa.x[array_ind];

      if (C_axis == zero_plane) {
        update_sign[0] = -1.;
      } else {
        update_sign[1] = -1.;
        update_sign[2] = -1.;
      }

      E_s_component = (C_axis) ? &inputs.E_s.zx : &inputs.E_s.yx;
      J_c_component = (C_axis) ? &J_c.zx : &J_c.yx;
      J_s_component = (C_axis) ? &J_s.zx : &J_s.yx;
      break;
    case AxialDirection::Y:
      // Dealing with Jsource. i = cell_b, k = cell_c
      c_constant = inputs.C.b.y[array_ind];

      s_index.j -= inputs.I0.index;
      s_index.k -= inputs.K0.index;
      source_value = inputs.Jsource[s_index];

      cell_to_update = {cell_b, inputs.J0.index, cell_c};
      if (!zero_plane) { cell_to_update.j = inputs.J1.index; }

      conductive_aux = inputs.rho_cond.y[array_ind];
      kappa = inputs.matched_layer.kappa.y[array_ind];

      // again, odd - the others have a symmetric whereas here we have an
      // asymmetry
      if (C_axis != zero_plane) {
        update_sign[0] = -1.;
      } else {
        update_sign[1] = -1.;
      }
      if (C_axis) { update_sign[2] = -1.; }

      E_s_component = (C_axis) ? &inputs.E_s.zy : &inputs.E_s.xy;
      J_c_component = (C_axis) ? &J_c.zy : &J_c.xy;
      J_s_component = (C_axis) ? &J_s.zy : &J_s.xy;
      break;
    case AxialDirection::Z:
      // Dealing with Ksource. i = cell_b, j = cell_c
      c_constant = inputs.C.b.z[array_ind];

      s_index.j -= inputs.I0.index;
      s_index.k -= inputs.J0.index;
      source_value = inputs.Ksource[s_index];

      cell_to_update = {cell_b, cell_c, inputs.K0.index};
      if (!zero_plane) { cell_to_update.k = inputs.K1.index; }

      conductive_aux = inputs.rho_cond.z[array_ind];
      kappa = inputs.matched_layer.kappa.z[array_ind];

      if (C_axis != zero_plane) {
        update_sign[1] = -1.;
      } else {
        update_sign[0] = -1.;
        update_sign[2] = -1.;
      }

      E_s_component = (C_axis) ? &inputs.E_s.yz : &inputs.E_s.xz;
      J_c_component = (C_axis) ? &J_c.yz : &J_c.xz;
      J_s_component = (C_axis) ? &J_s.yz : &J_s.xz;
      break;
  }

  // we can now compute the field update values
  double E_split_update =
          c_constant * real(common_amplitude * common_phase * source_value);

  // update the relevant split-field component
  (*E_s_component)[cell_to_update] += update_sign[0] * E_split_update;
  // update the current density in a conductive medium
  if (is_conductive) {
    (*J_c_component)[cell_to_update] +=
            update_sign[1] * conductive_aux * E_split_update;
  }
  // update the current density in a dispersive medium
  if (inputs.params.is_disp_ml) {
    (*J_s_component)[cell_to_update] += update_sign[2] * kappa * gamma /
                                        (2 * inputs.params.dt) * E_split_update;
  }
}
