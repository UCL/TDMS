#include "simulation_manager/simulation_manager.h"

#include "cell_coordinate.h"

using namespace std;
using namespace tdms_phys_constants;
using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;

void SimulationManager::update_source_terms_pulsed(
        double time_H, bool is_conductive, CurrentDensitySplitField &J_c,
        CurrentDensitySplitField &J_s) {
  //! Exit now if Ksource is empty, to avoid seg-faults
  if (inputs.Ksource.is_empty()) { return; }

  //! Simulation dimensions
  int I_tot = n_Yee_cells().i, J_tot = n_Yee_cells().j;
  //! The material constant that appears in the update equation
  double c_constant = inputs.C.b.z[inputs.K0.index];
  //! Conductive aux of the cell being updated
  double conductive_aux = inputs.rho_cond.z[inputs.K0.index];
  //! Dispersion prefactor in update equations
  double dispersion_factor = inputs.matched_layer.kappa.z[inputs.K0.index] *
                             inputs.matched_layer.gamma[inputs.K0.index] /
                             (2. * inputs.params.dt);
  //! Common amplitude factor in update equations
  double common_amplitude =
          c_constant * exp(-1.0 * DCPI *
                           pow((time_H - inputs.params.to_l +
                                inputs.params.delta.dz / LIGHT_V / 2.) /
                                       inputs.params.hwhm,
                               2));
  //! Common phase term in update equations
  complex<double> common_phase =
          -1.0 * IMAGINARY_UNIT *
          exp(-IMAGINARY_UNIT *
              fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                   2. * DCPI));
  //! Source term index to read from
  SourceIndex s_index;
  //! split-field cell index to update
  CellCoordinate cell_to_update;
  //! The update to apply to the split-field
  double split_field_update;

  if (J_tot == 0) {
    int j = 0;
    for (int i = 0; i < I_tot + 1; i++) {
      s_index = {2, i - inputs.I0.index, j};
      cell_to_update = {i, j, inputs.K0.index};

      split_field_update =
              common_amplitude * real(inputs.Ksource[s_index] * common_phase);

      inputs.E_s.yz[cell_to_update] -= split_field_update;
      if (is_conductive) {
        J_c.yz[cell_to_update] += conductive_aux * split_field_update;
      }
      if (inputs.params.is_disp_ml) {
        J_s.yz[cell_to_update] -= dispersion_factor * split_field_update;
      }
    }
  } else {
    for (int j = 0; j < J_tot; j++) {
      for (int i = 0; i < I_tot + 1; i++) {
        s_index = {2, i - inputs.I0.index, j - inputs.J0.index};
        cell_to_update = {i, j, inputs.K0.index};

        split_field_update =
                common_amplitude * real(inputs.Ksource[s_index] * common_phase);

        inputs.E_s.yz[cell_to_update] -= split_field_update;
        if (is_conductive) {
          J_c.yz[cell_to_update] += conductive_aux * split_field_update;
        }
        if (inputs.params.is_disp_ml) {
          J_s.yz[cell_to_update] -= dispersion_factor * split_field_update;
        }
      }
    }
  }

  for (int j = 0; j < J_tot + 1; j++) {
    for (int i = 0; i < I_tot; i++) {
      s_index = {3, i - inputs.I0.index, j - inputs.J0.index};
      cell_to_update = {i, j, inputs.K0.index};

      split_field_update =
              common_amplitude * real(inputs.Ksource[s_index] * common_phase);

      inputs.E_s.xz[cell_to_update] += split_field_update;
      if (is_conductive) {
        J_c.xz[cell_to_update] -= conductive_aux * split_field_update;
      }
      if (inputs.params.is_disp_ml) {
        J_s.xz[cell_to_update] += dispersion_factor * split_field_update;
      }
    }
  }
}

void SimulationManager::update_source_terms_pulsed(double time_E, int tind) {
  // Exit if there are no sources to update
  if (inputs.Ksource.is_empty()) { return; }

  //! Simulation dimensions
  int I_tot = n_Yee_cells().i, J_tot = n_Yee_cells().j;
  //! The material constant that appears in the update equation
  double d_constant = inputs.D.b.z[inputs.K0.index - 1];
  //! Common amplitude factor in update equations
  double common_amplitude =
          exp(-1. * DCPI *
              pow((time_E - inputs.params.to_l) / inputs.params.hwhm, 2.));
  //! Common phase term in update equations
  complex<double> common_phase =
          -1. * IMAGINARY_UNIT *
          exp(-1. * IMAGINARY_UNIT *
              fmod(inputs.params.omega_an * (time_E - inputs.params.to_l),
                   2 * DCPI));
  //! Source term index to read from
  SourceIndex s_index;
  //! split-field cell index to update
  CellCoordinate cell_to_update;

  if (J_tot == 0) {
    int j = 0;
    for (int i = 0; i < I_tot + 1; i++) {
      s_index = {1, i - inputs.I0.index, j};
      cell_to_update = {i, j, inputs.K0.index - 1};

      // Update magnetic split field
      inputs.H_s.xz[cell_to_update] -=
              d_constant * common_amplitude *
              real(inputs.Ksource[s_index] * common_phase);
      // Update broadband source term
      if (inputs.params.eyi_present) {
        inputs.H_s.xz[cell_to_update] -= d_constant * inputs.Ei.y[tind][j][i];
      }
    }
    for (int i = 0; i < I_tot; i++) {
      s_index = {0, i - inputs.I0.index, j};
      cell_to_update = {i, j, inputs.K0.index - 1};

      // Update magnetic split field
      inputs.H_s.yz[cell_to_update] +=
              d_constant * common_amplitude *
              real(inputs.Ksource[s_index] * common_phase);
      // Update broadband source term
      if (inputs.params.exi_present) {
        inputs.H_s.yz[cell_to_update] += d_constant * inputs.Ei.x[tind][j][i];
      }
    }
  } else {
    for (int j = 0; j < J_tot; j++) {
      for (int i = 0; i < I_tot + 1; i++) {
        s_index = {1, i - inputs.I0.index, j - inputs.J0.index};
        cell_to_update = {i, j, inputs.K0.index - 1};

        inputs.H_s.xz[cell_to_update] -=
                d_constant * common_amplitude *
                real(inputs.Ksource[s_index] * common_phase);
        if (inputs.params.eyi_present) {
          inputs.H_s.xz[cell_to_update] -= d_constant * inputs.Ei.y[tind][j][i];
        }
      }
    }
    for (int j = 0; j < J_tot + 1; j++) {
      for (int i = 0; i < I_tot; i++) {
        s_index = {0, i - inputs.I0.index, j - inputs.J0.index};
        cell_to_update = {i, j, inputs.K0.index};

        inputs.H_s.yz[cell_to_update] +=
                d_constant * common_amplitude *
                real(inputs.Ksource[s_index] * common_phase);
        if (inputs.params.exi_present) {
          inputs.H_s.yz[cell_to_update] += d_constant * inputs.Ei.x[tind][j][i];
        }
      }
    }
  }
}
