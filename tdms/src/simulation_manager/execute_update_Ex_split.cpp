#include "simulation_manager/simulation_manager.h"

#include <omp.h>

#include "numerical_derivative.h"

using namespace tdms_phys_constants;
using namespace std;

void SimulationManager::update_Exy(LoopVariables &lv) {
  // Fetch simulation dimensions for easy access
  int I_tot = n_Yee_cells().i, J_tot = n_Yee_cells().j, K_tot = n_Yee_cells().k;
  // Get the thread number
  int n = omp_get_thread_num();

#pragma omp for
  for (int k = 0; k < (K_tot + 1); k++) {
    for (int i = 0; i < I_tot; i++) {
      for (int j = 1; j < J_tot; j++) {
        double rho = 0.;
        int k_loc = k, array_ind = 0;
        double Ca, Cb, Cc;

        if (inputs.params.is_structure) {
          if (k > inputs.params.pml.Dzl &&
              k < (inputs.params.pml.Dzl + lv.n_non_pml_cells_in_K)) {
            if ((k - inputs.structure[i][1]) <
                        (lv.n_non_pml_cells_in_K + inputs.params.pml.Dzl) &&
                (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
              k_loc = k - inputs.structure[i][1];
            else if ((k - inputs.structure[i][1]) >=
                     (lv.n_non_pml_cells_in_K + inputs.params.pml.Dzl))
              k_loc = inputs.params.pml.Dzl + lv.n_non_pml_cells_in_K - 1;
            else
              k_loc = inputs.params.pml.Dzl + 1;
          }
        }
        if (!inputs.params.is_multilayer) {
          array_ind = j;
        } else {
          array_ind = (J_tot + 1) * k_loc + j;
        }

        // use the average of material parameters between nodes
        if (inputs.materials[k][j][i] || inputs.materials[k][j][i + 1]) {
          rho = 0.;
          if (!inputs.materials[k][j][i]) {
            Ca = inputs.C.a.y[array_ind];
            Cb = inputs.C.b.y[array_ind];
            if (inputs.params.is_disp_ml) Cc = inputs.C.c.y[array_ind];
            else
              Cc = 0.;
          } else {
            Ca = inputs.Cmaterial.a.y[inputs.materials[k][j][i] - 1];
            Cb = inputs.Cmaterial.b.y[inputs.materials[k][j][i] - 1];
            Cc = inputs.Cmaterial.c.y[inputs.materials[k][j][i] - 1];
          }

          if (inputs.params.interp_mat_props) {
            if (!inputs.materials[k][j][i + 1]) {
              Ca = Ca + inputs.C.a.y[array_ind];
              Cb = Cb + inputs.C.b.y[array_ind];
              if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.y[array_ind];
            } else {
              Ca = Ca + inputs.Cmaterial.a.y[inputs.materials[k][j][i + 1] - 1];
              Cb = Cb + inputs.Cmaterial.b.y[inputs.materials[k][j][i + 1] - 1];
              Cc = Cc + inputs.Cmaterial.c.y[inputs.materials[k][j][i + 1] - 1];
            }
            Ca = Ca / 2.;
            Cb = Cb / 2.;
            Cc = Cc / 2.;
          }
        } else {
          Ca = inputs.C.a.y[array_ind];
          Cb = inputs.C.b.y[array_ind];
          if (inputs.params.is_disp_ml) Cc = inputs.C.c.y[array_ind];
          else
            Cc = 0.;
          if (lv.is_conductive) rho = inputs.rho_cond.y[array_ind];
        }

        double alpha_l = 0.;
        double beta_l = 0.;
        double gamma_l = 0.;
        double kappa_l = 1.;
        double sigma_l = 0.;

        if (lv.is_dispersive || inputs.params.is_disp_ml) {
          sigma_l = inputs.matched_layer.sigma.y[array_ind];
          kappa_l = inputs.matched_layer.kappa.y[array_ind];
          alpha_l = inputs.matched_layer.alpha[k_loc];
          beta_l = inputs.matched_layer.beta[k_loc];
          gamma_l = inputs.matched_layer.gamma[k_loc];
          if (inputs.materials[k][j][i] || inputs.materials[k][j][i + 1]) {
            if (inputs.materials[k][j][i]) {
              alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
              beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
              gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
            } else {
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
            }

            if (inputs.materials[k][j][i + 1]) {
              alpha_l += inputs.alpha[inputs.materials[k][j][i + 1] - 1];
              beta_l += inputs.beta[inputs.materials[k][j][i + 1] - 1];
              gamma_l += inputs.gamma[inputs.materials[k][j][i + 1] - 1];
            } else {
              alpha_l += inputs.matched_layer.alpha[k_loc];
              beta_l += inputs.matched_layer.beta[k_loc];
              gamma_l += inputs.matched_layer.gamma[k_loc];
            }
            alpha_l = alpha_l / 2.;
            beta_l = beta_l / 2.;
            gamma_l = gamma_l / 2.;
          }
        }

        double Enp1, Jnp1;
        if (solver_method == SolverMethod::FiniteDifference) {
          Enp1 = Ca * inputs.E_s.xy[k][j][i] +
                 Cb * (inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i] -
                       inputs.H_s.zy[k][j - 1][i] - inputs.H_s.zx[k][j - 1][i]);
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.xy[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dy *
                            ((1 + alpha_l) * lv.J_s.xy[k][j][i] +
                             beta_l * lv.J_nm1.xy[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dy * lv.J_c.xy[k][j][i];
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l) {
            Jnp1 = alpha_l * lv.J_s.xy[k][j][i] +
                   beta_l * lv.J_nm1.xy[k][j][i] +
                   kappa_l * gamma_l / (2. * inputs.params.dt) *
                           (Enp1 - lv.E_nm1.xy[k][j][i]);
            Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xy[k][j][i];

            lv.E_nm1.xy[k][j][i] = inputs.E_s.xy[k][j][i];
            lv.J_nm1.xy[k][j][i] = lv.J_s.xy[k][j][i];
            lv.J_s.xy[k][j][i] = Jnp1;
          }

          if (lv.is_conductive && rho) {
            lv.J_c.xy[k][j][i] -= rho * (Enp1 + inputs.E_s.xy[k][j][i]);
          }

          inputs.E_s.xy[k][j][i] = Enp1;
        } else {// pseudo-spectral
          Enp1 = 0.0;
          // Enp1 = Ca*E_s.xy[k][j][i]+Cb*(H_s.zy[k][j][i] + H_s.zx[k][j][i] -
          // H_s.zy[k][j-1][i] - H_s.zx[k][j-1][i]);
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.xy[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dy *
                            ((1 + alpha_l) * lv.J_s.xy[k][j][i] +
                             beta_l * lv.J_nm1.xy[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dy * lv.J_c.xy[k][j][i];
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l) {
            Jnp1 = alpha_l * lv.J_s.xy[k][j][i] +
                   beta_l * lv.J_nm1.xy[k][j][i] +
                   kappa_l * gamma_l / (2. * inputs.params.dt) *
                           (Enp1 - lv.E_nm1.xy[k][j][i]);
            Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xy[k][j][i];

            lv.E_nm1.xy[k][j][i] = inputs.E_s.xy[k][j][i];
            lv.J_nm1.xy[k][j][i] = lv.J_s.xy[k][j][i];
            lv.J_s.xy[k][j][i] = Jnp1;
          }

          if (lv.is_conductive && rho) {
            lv.J_c.xy[k][j][i] -= rho * (Enp1 + inputs.E_s.xy[k][j][i]);
          }

          eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
          eh_vec[n][j][1] = 0.;
          PSTD.ca[n][j - 1] = Ca;
          PSTD.cb[n][j - 1] = Cb;
        }
      }
      if (solver_method == SolverMethod::PseudoSpectral && J_tot > 1) {
        int j = 0;
        eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
        eh_vec[n][j][1] = 0.;
        first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ey, PSTD.N_ey,
                         inputs.E_s.xy.plan_f[n], inputs.E_s.xy.plan_b[n]);
        for (j = 1; j < J_tot; j++) {
          inputs.E_s.xy[k][j][i] =
                  PSTD.ca[n][j - 1] * inputs.E_s.xy[k][j][i] +
                  PSTD.cb[n][j - 1] * eh_vec[n][j][0] / ((double) PSTD.N_ey);
        }
      }
    }
  }
}

void SimulationManager::update_Exz(LoopVariables &lv) {
  // Fetch simulation dimensions for easy access
  int I_tot = n_Yee_cells().i, K_tot = n_Yee_cells().k;
  // Get the thread number
  int n = omp_get_thread_num();

#pragma omp for
  for (int j = 0; j < lv.J_loop_upper_bound_plus_1; j++) {
    for (int i = 0; i < I_tot; i++) {
      for (int k = 1; k < K_tot; k++) {
        double rho = 0.;
        int k_loc = k;
        double Ca, Cb, Cc;

        if (inputs.params.is_structure)
          if (k > inputs.params.pml.Dzl &&
              k < (inputs.params.pml.Dzl + lv.n_non_pml_cells_in_K)) {
            if ((k - inputs.structure[i][1]) <
                        (lv.n_non_pml_cells_in_K + inputs.params.pml.Dzl) &&
                (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
              k_loc = k - inputs.structure[i][1];
            else if ((k - inputs.structure[i][1]) >=
                     (lv.n_non_pml_cells_in_K + inputs.params.pml.Dzl))
              k_loc = inputs.params.pml.Dzl + lv.n_non_pml_cells_in_K - 1;
            else
              k_loc = inputs.params.pml.Dzl + 1;
          }
        // use the average of material parameters between nodes
        if (inputs.materials[k][j][i] || inputs.materials[k][j][i + 1]) {
          rho = 0.;
          if (!inputs.materials[k][j][i]) {
            Ca = inputs.C.a.z[k_loc];
            Cb = inputs.C.b.z[k_loc];
            if (inputs.params.is_disp_ml) Cc = inputs.C.c.z[k_loc];
            else
              Cc = 0.;
          } else {
            Ca = inputs.Cmaterial.a.z[inputs.materials[k][j][i] - 1];
            Cb = inputs.Cmaterial.b.z[inputs.materials[k][j][i] - 1];
            Cc = inputs.Cmaterial.c.z[inputs.materials[k][j][i] - 1];
          }

          if (inputs.params.interp_mat_props) {
            if (!inputs.materials[k][j][i + 1]) {
              Ca = Ca + inputs.C.a.z[k_loc];
              Cb = Cb + inputs.C.b.z[k_loc];
              if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.z[k_loc];
            } else {
              Ca = Ca + inputs.Cmaterial.a.z[inputs.materials[k][j][i + 1] - 1];
              Cb = Cb + inputs.Cmaterial.b.z[inputs.materials[k][j][i + 1] - 1];
              Cc = Cc + inputs.Cmaterial.c.z[inputs.materials[k][j][i + 1] - 1];
            }
            Ca = Ca / 2.;
            Cb = Cb / 2.;
            Cc = Cc / 2.;
          }
        } else {
          Ca = inputs.C.a.z[k_loc];
          Cb = inputs.C.b.z[k_loc];
          if (inputs.params.is_disp_ml) Cc = inputs.C.c.z[k_loc];
          else
            Cc = 0.;
          if (lv.is_conductive) rho = inputs.rho_cond.z[k_loc];
        }

        double alpha_l = 0.;
        double beta_l = 0.;
        double gamma_l = 0.;
        double kappa_l = 1.;
        double sigma_l = 0.;

        if (lv.is_dispersive || inputs.params.is_disp_ml) {
          sigma_l = inputs.matched_layer.sigma.z[k_loc];
          kappa_l = inputs.matched_layer.kappa.z[k_loc];
          alpha_l = inputs.matched_layer.alpha[k_loc];
          beta_l = inputs.matched_layer.beta[k_loc];
          gamma_l = inputs.matched_layer.gamma[k_loc];
          if (inputs.materials[k][j][i] || inputs.materials[k][j][i + 1]) {
            if (inputs.materials[k][j][i]) {
              alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
              beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
              gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
            } else {
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
            }

            if (inputs.materials[k][j][i + 1]) {
              alpha_l += inputs.alpha[inputs.materials[k][j][i + 1] - 1];
              beta_l += inputs.beta[inputs.materials[k][j][i + 1] - 1];
              gamma_l += inputs.gamma[inputs.materials[k][j][i + 1] - 1];
            } else {
              alpha_l += inputs.matched_layer.alpha[k_loc];
              beta_l += inputs.matched_layer.beta[k_loc];
              gamma_l += inputs.matched_layer.gamma[k_loc];
            }
            alpha_l = alpha_l / 2.;
            beta_l = beta_l / 2.;
            gamma_l = gamma_l / 2.;
          }
        }
        // only things used here are the greek letters and the Ca, Cb, Cc
        // variables...
        double Enp1, Jnp1;
        if (solver_method == SolverMethod::FiniteDifference) {
          Enp1 = Ca * inputs.E_s.xz[k][j][i] +
                 Cb * (inputs.H_s.yx[k - 1][j][i] + inputs.H_s.yz[k - 1][j][i] -
                       inputs.H_s.yx[k][j][i] - inputs.H_s.yz[k][j][i]);
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.xz[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dz *
                            ((1 + alpha_l) * lv.J_s.xz[k][j][i] +
                             beta_l * lv.J_nm1.xz[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dz * lv.J_c.xz[k][j][i];
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l) {
            Jnp1 = alpha_l * lv.J_s.xz[k][j][i] +
                   beta_l * lv.J_nm1.xz[k][j][i] +
                   kappa_l * gamma_l / (2. * inputs.params.dt) *
                           (Enp1 - lv.E_nm1.xz[k][j][i]);
            Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xz[k][j][i];
            lv.E_nm1.xz[k][j][i] = inputs.E_s.xz[k][j][i];
            lv.J_nm1.xz[k][j][i] = lv.J_s.xz[k][j][i];
            lv.J_s.xz[k][j][i] = Jnp1;
          }

          if (lv.is_conductive && rho) {
            lv.J_c.xz[k][j][i] -= rho * (Enp1 + inputs.E_s.xz[k][j][i]);
          }

          inputs.E_s.xz[k][j][i] = Enp1;
        } else {// psuedo-spectral
          // Enp1 = Ca*E_s.xz[k][j][i]+Cb*(H_s.yx[k-1][j][i] + H_s.yz[k-1][j][i]
          // - H_s.yx[k][j][i] - H_s.yz[k][j][i]);
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.xz[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dz *
                            ((1 + alpha_l) * lv.J_s.xz[k][j][i] +
                             beta_l * lv.J_nm1.xz[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dz * lv.J_c.xz[k][j][i];
          if ((lv.is_dispersive || inputs.params.is_disp_ml) && gamma_l) {
            Jnp1 = alpha_l * lv.J_s.xz[k][j][i] +
                   beta_l * lv.J_nm1.xz[k][j][i] +
                   kappa_l * gamma_l / (2. * inputs.params.dt) *
                           (Enp1 - lv.E_nm1.xz[k][j][i]);
            Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xz[k][j][i];
            lv.E_nm1.xz[k][j][i] = inputs.E_s.xz[k][j][i];
            lv.J_nm1.xz[k][j][i] = lv.J_s.xz[k][j][i];
            lv.J_s.xz[k][j][i] = Jnp1;
          }

          if (lv.is_conductive && rho) {
            lv.J_c.xz[k][j][i] -= rho * (Enp1 + inputs.E_s.xz[k][j][i]);
          }

          eh_vec[n][k][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
          eh_vec[n][k][1] = 0.;
          PSTD.ca[n][k - 1] = Ca;
          PSTD.cb[n][k - 1] = Cb;
        }
      }
      if (solver_method == SolverMethod::PseudoSpectral) {
        int k = 0;
        eh_vec[n][k][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
        eh_vec[n][k][1] = 0.;

        first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ez, PSTD.N_ez,
                         inputs.E_s.xz.plan_f[n], inputs.E_s.xz.plan_b[n]);

        for (k = 1; k < K_tot; k++) {
          inputs.E_s.xz[k][j][i] =
                  PSTD.ca[n][k - 1] * inputs.E_s.xz[k][j][i] -
                  PSTD.cb[n][k - 1] * eh_vec[n][k][0] / ((double) PSTD.N_ez);
        }
      }
    }
  }
}
