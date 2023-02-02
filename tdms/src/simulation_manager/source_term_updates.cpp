#include "simulation_manager/simulation_manager.h"

#include "cell_coordinate.h"

using namespace std;
using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;

void SimulationManager::update_Isource_terms_steadystate(
        double time_H, bool is_conductive, CurrentDensitySplitField &J_c,
        CurrentDensitySplitField &J_s) {
  // Only run update equations is source data was provided
  if (inputs.Isource.is_empty()) { return; }

  int array_ind;
  int I_tot = n_Yee_cells().i;

  complex<double> commonPhase = exp(
          -IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_H, 2. * DCPI));
  double commonAmplitude = linear_ramp(time_H);

  // Update across I0, provided a source term is here
  if (inputs.I0.apply) {
    for (int k = (inputs.K0.index); k <= (inputs.K1.index); k++) {
      for (int j = (inputs.J0.index); j <= (inputs.J1.index); j++) {
        if (!inputs.params.is_multilayer) {
          array_ind = inputs.I0.index;
        } else {
          array_ind = (I_tot + 1) * k + inputs.I0.index;
        }

        // Update zx
        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          // Compute common prefactor
          SourceIndex s_index{2, j - inputs.J0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.x[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Isource[s_index]);

          inputs.E_s.zx[k][j][inputs.I0.index] -= common_factor;
          if (is_conductive) {
            J_c.zx[k][j][inputs.I0.index] +=
                    inputs.rho_cond.x[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.zx[k][j][inputs.I0.index] +=
                    inputs.matched_layer.kappa.x[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }

        // Update yx
        if (j < (inputs.J1.index)) {
          // Compute common prefactor
          SourceIndex s_index{3, j - inputs.J0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.x[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Isource[s_index]);

          inputs.E_s.yx[k][j][inputs.I0.index] += common_factor;
          if (is_conductive) {
            J_c.yx[k][j][inputs.I0.index] -=
                    inputs.rho_cond.x[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.yx[k][j][inputs.I0.index] -=
                    inputs.matched_layer.kappa.x[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
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
          array_ind = (I_tot + 1) * k + inputs.I1.index;
        }

        // Update zx
        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          // Common prefactors
          SourceIndex s_index{6, j - inputs.J0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.x[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Isource[s_index]);

          inputs.E_s.zx[k][j][inputs.I1.index] += common_factor;
          if (is_conductive) {
            J_c.zx[k][j][inputs.I1.index] -=
                    inputs.rho_cond.x[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.zx[k][j][inputs.I1.index] -=
                    inputs.matched_layer.kappa.x[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }

        if (j < (inputs.J1.index)) {
          // Common prefactors
          SourceIndex s_index{7, j - inputs.J0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.x[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Isource[s_index]);

          inputs.E_s.yx[k][j][inputs.I1.index] -= common_factor;
          if (is_conductive) {
            J_c.yx[k][j][inputs.I1.index] +=
                    inputs.rho_cond.x[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.yx[k][j][inputs.I1.index] +=
                    inputs.matched_layer.kappa.x[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
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
  int J_tot = n_Yee_cells().j;

  complex<double> commonPhase = exp(
          -IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_H, 2. * DCPI));
  double commonAmplitude = linear_ramp(time_H);

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
            array_ind = (J_tot + 1) * k + inputs.J0.index;
          }

          SourceIndex s_index{2, i - inputs.I0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.y[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Jsource[s_index]);

          inputs.E_s.zy[k][inputs.J0.index][i] += common_factor;
          if (is_conductive) {
            J_c.zy[k][inputs.J0.index][i] -=
                    inputs.rho_cond.y[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.zy[k][inputs.J0.index][i] -=
                    inputs.matched_layer.kappa.y[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }

        if (i < inputs.I1.index) {
          SourceIndex s_index{3, i - inputs.I0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.y[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Jsource[s_index]);

          inputs.E_s.xy[k][inputs.J0.index][i] -= common_factor;
          if (is_conductive) {
            J_c.xy[k][inputs.J0.index][i] +=
                    inputs.rho_cond.y[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.xy[k][(inputs.J0.index)][i] +=
                    inputs.matched_layer.kappa.y[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
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
          array_ind = (J_tot + 1) * k + inputs.J1.index;
        }

        if (k < (inputs.K1.index) ||
            inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
          SourceIndex s_index{6, i - inputs.I0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.y[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Jsource[s_index]);

          inputs.E_s.zy[k][(inputs.J1.index)][i] -= common_factor;
          if (is_conductive) {
            J_c.zy[k][(inputs.J1.index)][i] +=
                    inputs.rho_cond.y[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.zy[k][(inputs.J1.index)][i] -=
                    inputs.matched_layer.kappa.y[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }

        if (i < (inputs.I1.index)) {
          SourceIndex s_index{7, i - inputs.I0.index, k - inputs.K0.index};
          double common_factor =
                  inputs.C.b.y[array_ind] *
                  real(commonAmplitude * commonPhase * inputs.Jsource[s_index]);

          inputs.E_s.xy[k][(inputs.J1.index)][i] += common_factor;
          if (is_conductive) {
            J_c.xy[k][(inputs.J1.index)][i] -=
                    inputs.rho_cond.y[array_ind] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.xy[k][(inputs.J1.index)][i] +=
                    inputs.matched_layer.kappa.y[array_ind] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
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

  // TODO: Line of doom hangover! k is effectively this at all times, due to
  // previous scoping errors
  int k = inputs.K1.index;

  complex<double> commonPhase = exp(
          -IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_H, 2. * DCPI));
  double commonAmplitude = linear_ramp(time_H);

  // Update across K0, provided a source term is there
  if (inputs.K0.apply) {
    for (int j = inputs.J0.index; j <= inputs.J1.index; j++) {
      for (int i = inputs.I0.index; i <= inputs.I1.index; i++) {
        if (j < (inputs.J1.index)) {
          SourceIndex s_index{2, i - inputs.I0.index, j - inputs.J0.index};
          double common_factor =
                  inputs.C.b.z[inputs.K0.index] *
                  real(commonAmplitude * commonPhase * inputs.Ksource[s_index]);

          inputs.E_s.yz[(inputs.K0.index)][j][i] -= common_factor;
          if (is_conductive) {
            J_c.yz[(inputs.K0.index)][j][i] +=
                    inputs.rho_cond.z[(inputs.K0.index)] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.yz[(inputs.K0.index)][j][i] -=
                    inputs.matched_layer.kappa.z[(inputs.K0.index)] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }

        if (i < (inputs.I1.index)) {
          SourceIndex s_index{3, i - inputs.I0.index, j - inputs.J0.index};
          double common_factor =
                  inputs.C.b.z[inputs.K0.index] *
                  real(commonAmplitude * commonPhase * inputs.Ksource[s_index]);

          inputs.E_s.xz[(inputs.K0.index)][j][i] += common_factor;
          if (is_conductive) {
            J_c.xz[(inputs.K0.index)][j][i] -=
                    inputs.rho_cond.z[(inputs.K0.index)] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.xz[(inputs.K0.index)][j][i] +=
                    inputs.matched_layer.kappa.z[(inputs.K0.index)] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }
      }
    }
  }

  // Update across K1, provided a source term is there
  if (inputs.K1.apply) {
    for (int j = inputs.J0.index; j <= inputs.J1.index; j++) {
      for (int i = inputs.I0.index; i <= inputs.I1.index; i++) {
        if (j < (inputs.J1.index)) {
          SourceIndex s_index{6, i - inputs.I0.index, j - inputs.J0.index};
          double common_factor =
                  inputs.C.b.z[inputs.K1.index] *
                  real(commonAmplitude * commonPhase * inputs.Ksource[s_index]);

          inputs.E_s.yz[(inputs.K1.index)][j][i] += common_factor;
          if (is_conductive) {
            J_c.yz[inputs.K1.index][j][i] -=
                    inputs.rho_cond.z[(inputs.K1.index)] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.yz[(inputs.K1.index)][j][i] +=
                    inputs.matched_layer.kappa.z[(inputs.K1.index)] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }

        if (i < (inputs.I1.index)) {
          SourceIndex s_index{7, i - inputs.I0.index, j - inputs.J0.index};
          double common_factor =
                  inputs.C.b.z[inputs.K1.index] *
                  real(commonAmplitude * commonPhase * inputs.Ksource[s_index]);

          inputs.E_s.xz[(inputs.K1.index)][j][i] -= common_factor;
          if (is_conductive) {
            J_c.xz[(inputs.K1.index)][j][i] +=
                    inputs.rho_cond.z[(inputs.K1.index)] * common_factor;
          }
          if (inputs.params.is_disp_ml) {
            J_s.xz[(inputs.K1.index)][j][i] -=
                    inputs.matched_layer.kappa.z[(inputs.K1.index)] *
                    inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                    common_factor;
          }
        }
      }
    }
  }
}
