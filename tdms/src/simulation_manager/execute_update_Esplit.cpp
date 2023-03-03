#include "simulation_manager.h"

#include <omp.h>

#include "numerical_derivative.h"

using namespace tdms_phys_constants;
using namespace std;

void SimulationManager::iterate_E_split(LoopVariables &lv) {
  // Fetch simulation dimensions for easy asscess
  int I_tot = n_Yee_cells().I_tot(), J_tot = n_Yee_cells().J_tot(),
      K_tot = n_Yee_cells().K_tot();
  // Get the thread number
  int n = omp_get_thread_num();

  // remove the loop variables once the refactor is done...
  int i, j, k;

  double rho;
  double alpha_l, beta_l,
          gamma_l;//< alpha, beta, gamma parameters of the layer the local
                  // thread is examining
  double kappa_l, sigma_l;//< kappa, sigma parameters of the layer the local
                          // thread is examining
  int k_loc, array_ind;
  double Ca, Cb, Cc;// used by interpolation scheme
                    // the C and D vars for free space and pml
  double Enp1, Jnp1;

  if (inputs.params.dimension == THREE ||
      inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
    if (solver_method == SolverMethod::FiniteDifference) {
      // FDTD, E_s.xy
#pragma omp for
      for (k = 0; k < (K_tot + 1); k++)
        for (j = 1; j < J_tot; j++)
          for (i = 0; i < I_tot; i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (!inputs.params.is_multilayer) array_ind = j;
            else
              array_ind = (J_tot + 1) * k_loc + j;

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
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.y[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a.y[inputs.materials[k][j][i + 1] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.y[inputs.materials[k][j][i + 1] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.y[inputs.materials[k][j][i + 1] - 1];
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

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
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


            Enp1 = Ca * inputs.E_s.xy[k][j][i] +
                   Cb * (inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i] -
                         inputs.H_s.zy[k][j - 1][i] -
                         inputs.H_s.zx[k][j - 1][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.xy[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dy *
                              ((1 + alpha_l) * lv.J_s.xy[k][j][i] +
                               beta_l * lv.J_nm1.xy[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dy * lv.J_c.xy[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
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
          }
      // FDTD, E_s.xy
    } else {
#pragma omp for
      for (k = 0; k < (K_tot + 1); k++)
        for (i = 0; i < I_tot; i++) {
          for (j = 1; j < J_tot; j++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (!inputs.params.is_multilayer) array_ind = j;
            else
              array_ind = (J_tot + 1) * k_loc + j;

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
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.y[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a.y[inputs.materials[k][j][i + 1] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.y[inputs.materials[k][j][i + 1] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.y[inputs.materials[k][j][i + 1] - 1];
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

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
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


            Enp1 = 0.0;
            // Enp1 = Ca*E_s.xy[k][j][i]+Cb*(H_s.zy[k][j][i] + H_s.zx[k][j][i] -
            // H_s.zy[k][j-1][i] - H_s.zx[k][j-1][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.xy[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dy *
                              ((1 + alpha_l) * lv.J_s.xy[k][j][i] +
                               beta_l * lv.J_nm1.xy[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dy * lv.J_c.xy[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
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
          if (J_tot > 1) {
            j = 0;
            eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
            eh_vec[n][j][1] = 0.;
            first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ey, PSTD.N_ey,
                             inputs.E_s.xy.plan_f[n], inputs.E_s.xy.plan_b[n]);
            for (j = 1; j < J_tot; j++) {
              inputs.E_s.xy[k][j][i] =
                      PSTD.ca[n][j - 1] * inputs.E_s.xy[k][j][i] +
                      PSTD.cb[n][j - 1] * eh_vec[n][j][0] /
                              ((double) PSTD.N_ey);
            }
          }
        }
      // PSTD, E_s.xy
    }// if (solver_method == DerivativeMethod::FiniteDifference) (else
     // PseudoSpectral)
    // E_s.xz updates
    if (solver_method == SolverMethod::FiniteDifference) {
#pragma omp for
      for (k = 1; k < K_tot; k++)
        for (j = 0; j < lv.J_tot_p1_bound; j++)
          for (i = 0; i < I_tot; i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
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
                  Ca = Ca +
                       inputs.Cmaterial.a.z[inputs.materials[k][j][i + 1] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.z[inputs.materials[k][j][i + 1] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.z[inputs.materials[k][j][i + 1] - 1];
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

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
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
            Enp1 = Ca * inputs.E_s.xz[k][j][i] +
                   Cb * (inputs.H_s.yx[k - 1][j][i] +
                         inputs.H_s.yz[k - 1][j][i] - inputs.H_s.yx[k][j][i] -
                         inputs.H_s.yz[k][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.xz[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dz *
                              ((1 + alpha_l) * lv.J_s.xz[k][j][i] +
                               beta_l * lv.J_nm1.xz[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dz * lv.J_c.xz[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
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
          }
      // FDTD, E_s.xz
    } else {
      //#pragma omp for
      for (j = 0; j < lv.J_tot_p1_bound; j++)
#pragma omp for
        for (i = 0; i < I_tot; i++) {
          for (k = 1; k < K_tot; k++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
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
                  Ca = Ca +
                       inputs.Cmaterial.a.z[inputs.materials[k][j][i + 1] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.z[inputs.materials[k][j][i + 1] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.z[inputs.materials[k][j][i + 1] - 1];
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

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
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
            // Enp1 = Ca*E_s.xz[k][j][i]+Cb*(H_s.yx[k-1][j][i] +
            // H_s.yz[k-1][j][i] - H_s.yx[k][j][i] - H_s.yz[k][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.xz[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dz *
                              ((1 + alpha_l) * lv.J_s.xz[k][j][i] +
                               beta_l * lv.J_nm1.xz[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dz * lv.J_c.xz[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
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
          k = 0;
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
      // PSTD, E_s.xz
    }// if (solver_method == DerivativeMethod::FiniteDifference) (else
     // PseudoSpectral)
    // E_s.yx updates
    if (solver_method == SolverMethod::FiniteDifference) {
      // FDTD, E_s.yx
#pragma omp for
      for (k = 0; k < (K_tot + 1); k++)
        for (j = 0; j < lv.J_tot_bound; j++)
          for (i = 1; i < I_tot; i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure) {
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            }
            if (!inputs.params.is_multilayer) array_ind = i;
            else
              array_ind = (I_tot + 1) * k_loc + i;

            // use the average of material parameters between nodes
            if (inputs.materials[k][j][i] ||
                inputs.materials[k][min(J_tot, j + 1)][i]) {
              rho = 0.;
              if (!inputs.materials[k][j][i]) {
                Ca = inputs.C.a.x[array_ind];
                Cb = inputs.C.b.x[array_ind];
                if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
                else
                  Cc = 0;
              } else {
                Ca = inputs.Cmaterial.a.x[inputs.materials[k][j][i] - 1];
                Cb = inputs.Cmaterial.b.x[inputs.materials[k][j][i] - 1];
                Cc = inputs.Cmaterial.c.x[inputs.materials[k][j][i] - 1];
              }
              if (inputs.params.interp_mat_props) {
                if (!inputs.materials[k][min(J_tot, j + 1)][i]) {
                  Ca = Ca + inputs.C.a.x[array_ind];
                  Cb = Cb + inputs.C.b.x[array_ind];
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.x[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a
                               .x[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cb = Cb +
                       inputs.Cmaterial.b
                               .x[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cc = Cc +
                       inputs.Cmaterial.c
                               .x[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                }

                Ca = Ca / 2.;
                Cb = Cb / 2.;
                Cc = Cc / 2.;
              }
            } else {
              Ca = inputs.C.a.x[array_ind];
              Cb = inputs.C.b.x[array_ind];
              if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
              else
                Cc = 0.;
              if (lv.is_conductive) rho = inputs.rho_cond.x[array_ind];
            }

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.x[array_ind];
              kappa_l = inputs.matched_layer.kappa.x[array_ind];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] ||
                  inputs.materials[k][min(J_tot, j + 1)][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k][min(J_tot, j + 1)][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
                  beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)]
                                                        [i] -
                                        1];
                  gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
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


            Enp1 = Ca * inputs.E_s.yx[k][j][i] +
                   Cb * (inputs.H_s.zx[k][j][i - 1] +
                         inputs.H_s.zy[k][j][i - 1] - inputs.H_s.zx[k][j][i] -
                         inputs.H_s.zy[k][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.yx[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dx *
                              ((1 + alpha_l) * lv.J_s.yx[k][j][i] +
                               beta_l * lv.J_nm1.yx[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dx * lv.J_c.yx[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.yx[k][j][i] +
                     beta_l * lv.J_nm1.yx[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.yx[k][j][i]);
              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx[k][j][i];
              lv.E_nm1.yx[k][j][i] = inputs.E_s.yx[k][j][i];
              lv.J_nm1.yx[k][j][i] = lv.J_s.yx[k][j][i];
              lv.J_s.yx[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.yx[k][j][i] -= rho * (Enp1 + inputs.E_s.yx[k][j][i]);
            }

            inputs.E_s.yx[k][j][i] = Enp1;
          }
      // FDTD, E_s.yx
    } else {
#pragma omp for
      for (k = 0; k < (K_tot + 1); k++)
        for (j = 0; j < lv.J_tot_bound; j++) {
          for (i = 1; i < I_tot; i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure) {
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            }
            if (!inputs.params.is_multilayer) array_ind = i;
            else
              array_ind = (I_tot + 1) * k_loc + i;

            // use the average of material parameters between nodes
            if (inputs.materials[k][j][i] ||
                inputs.materials[k][min(J_tot, j + 1)][i]) {
              rho = 0.;
              if (!inputs.materials[k][j][i]) {
                Ca = inputs.C.a.x[array_ind];
                Cb = inputs.C.b.x[array_ind];
                if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
                else
                  Cc = 0;
              } else {
                Ca = inputs.Cmaterial.a.x[inputs.materials[k][j][i] - 1];
                Cb = inputs.Cmaterial.b.x[inputs.materials[k][j][i] - 1];
                Cc = inputs.Cmaterial.c.x[inputs.materials[k][j][i] - 1];
              }
              if (inputs.params.interp_mat_props) {
                if (!inputs.materials[k][min(J_tot, j + 1)][i]) {
                  Ca = Ca + inputs.C.a.x[array_ind];
                  Cb = Cb + inputs.C.b.x[array_ind];
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.x[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a
                               .x[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cb = Cb +
                       inputs.Cmaterial.b
                               .x[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cc = Cc +
                       inputs.Cmaterial.c
                               .x[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                }

                Ca = Ca / 2.;
                Cb = Cb / 2.;
                Cc = Cc / 2.;
              }
            } else {
              Ca = inputs.C.a.x[array_ind];
              Cb = inputs.C.b.x[array_ind];
              if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
              else
                Cc = 0.;
              if (lv.is_conductive) rho = inputs.rho_cond.x[array_ind];
            }

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.x[array_ind];
              kappa_l = inputs.matched_layer.kappa.x[array_ind];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] ||
                  inputs.materials[k][min(J_tot, j + 1)][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k][min(J_tot, j + 1)][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
                  beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)]
                                                        [i] -
                                        1];
                  gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
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


            // Enp1 = Ca*E_s.yx[k][j][i]+Cb*(H_s.zx[k][j][i-1] +
            // H_s.zy[k][j][i-1] - H_s.zx[k][j][i] - H_s.zy[k][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.yx[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dx *
                              ((1 + alpha_l) * lv.J_s.yx[k][j][i] +
                               beta_l * lv.J_nm1.yx[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dx * lv.J_c.yx[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.yx[k][j][i] +
                     beta_l * lv.J_nm1.yx[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.yx[k][j][i]);
              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx[k][j][i];
              lv.E_nm1.yx[k][j][i] = inputs.E_s.yx[k][j][i];
              lv.J_nm1.yx[k][j][i] = lv.J_s.yx[k][j][i];
              lv.J_s.yx[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.yx[k][j][i] -= rho * (Enp1 + inputs.E_s.yx[k][j][i]);
            }

            eh_vec[n][i][0] = inputs.H_s.zx[k][j][i] + inputs.H_s.zy[k][j][i];
            eh_vec[n][i][1] = 0.;
            PSTD.ca[n][i - 1] = Ca;
            PSTD.cb[n][i - 1] = Cb;
          }
          i = 0;
          eh_vec[n][i][0] = inputs.H_s.zx[k][j][i] + inputs.H_s.zy[k][j][i];
          eh_vec[n][i][1] = 0.;

          first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ex, PSTD.N_ex,
                           inputs.E_s.yx.plan_f[n], inputs.E_s.yx.plan_b[n]);

          for (i = 1; i < I_tot; i++) {
            inputs.E_s.yx[k][j][i] =
                    PSTD.ca[n][i - 1] * inputs.E_s.yx[k][j][i] -
                    PSTD.cb[n][i - 1] * eh_vec[n][i][0] / ((double) PSTD.N_ex);
            // E_s.yx[k][j][i] = Enp1;
          }
        }
      // PSTD, E_s.yx
    }// if (solver_method == DerivativeMethod::FiniteDifference) (else
     // PseudoSpectral)
    // E_s.yz updates
    if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, E_s.yz
#pragma omp for
      for (k = 1; k < K_tot; k++)
        for (j = 0; j < lv.J_tot_bound; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (inputs.materials[k][j][i] ||
                inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                if (!inputs.materials[k][min(J_tot, j + 1)][i]) {
                  Ca = Ca + inputs.C.a.z[k_loc];
                  Cb = Cb + inputs.C.b.z[k_loc];
                  if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.z[k_loc];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a
                               .z[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cb = Cb +
                       inputs.Cmaterial.b
                               .z[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cc = Cc +
                       inputs.Cmaterial.c
                               .z[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
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

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.z[k_loc];
              kappa_l = inputs.matched_layer.kappa.z[k_loc];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] ||
                  inputs.materials[k][min(J_tot, j + 1)][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k][min(J_tot, j + 1)][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
                  beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)]
                                                        [i] -
                                        1];
                  gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
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
            Enp1 = Ca * inputs.E_s.yz[k][j][i] +
                   Cb * (inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i] -
                         inputs.H_s.xy[k - 1][j][i] -
                         inputs.H_s.xz[k - 1][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.yz[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dz *
                              ((1 + alpha_l) * lv.J_s.yz[k][j][i] +
                               beta_l * lv.J_nm1.yz[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dz * lv.J_c.yz[k][j][i];

            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.yz[k][j][i] +
                     beta_l * lv.J_nm1.yz[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.yz[k][j][i]);
              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz[k][j][i];
              lv.E_nm1.yz[k][j][i] = inputs.E_s.yz[k][j][i];
              lv.J_nm1.yz[k][j][i] = lv.J_s.yz[k][j][i];
              lv.J_s.yz[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.yz[k][j][i] -= rho * (Enp1 + inputs.E_s.yz[k][j][i]);
            }

            inputs.E_s.yz[k][j][i] = Enp1;
          }
      // FDTD, E_s.yz
    } else {
#pragma omp for
      for (j = 0; j < lv.J_tot_bound; j++)
        for (i = 0; i < (I_tot + 1); i++) {
          for (k = 1; k < K_tot; k++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (inputs.materials[k][j][i] ||
                inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                if (!inputs.materials[k][min(J_tot, j + 1)][i]) {
                  Ca = Ca + inputs.C.a.z[k_loc];
                  Cb = Cb + inputs.C.b.z[k_loc];
                  if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.z[k_loc];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a
                               .z[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cb = Cb +
                       inputs.Cmaterial.b
                               .z[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
                  Cc = Cc +
                       inputs.Cmaterial.c
                               .z[inputs.materials[k][min(J_tot, j + 1)][i] -
                                  1];
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

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.z[k_loc];
              kappa_l = inputs.matched_layer.kappa.z[k_loc];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] ||
                  inputs.materials[k][min(J_tot, j + 1)][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k][min(J_tot, j + 1)][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
                  beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)]
                                                        [i] -
                                        1];
                  gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)]
                                                          [i] -
                                          1];
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
            // Enp1 = Ca*E_s.yz[k][j][i]+Cb*(H_s.xy[k][j][i] + H_s.xz[k][j][i] -
            // H_s.xy[k-1][j][i] - H_s.xz[k-1][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.yz[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dz *
                              ((1 + alpha_l) * lv.J_s.yz[k][j][i] +
                               beta_l * lv.J_nm1.yz[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dz * lv.J_c.yz[k][j][i];

            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.yz[k][j][i] +
                     beta_l * lv.J_nm1.yz[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.yz[k][j][i]);
              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz[k][j][i];
              lv.E_nm1.yz[k][j][i] = inputs.E_s.yz[k][j][i];
              lv.J_nm1.yz[k][j][i] = lv.J_s.yz[k][j][i];
              lv.J_s.yz[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.yz[k][j][i] -= rho * (Enp1 + inputs.E_s.yz[k][j][i]);
            }

            eh_vec[n][k][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
            eh_vec[n][k][1] = 0.;
            PSTD.ca[n][k - 1] = Ca;
            PSTD.cb[n][k - 1] = Cb;
          }
          k = 0;
          eh_vec[n][k][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
          eh_vec[n][k][1] = 0.;
          first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ez, PSTD.N_ez,
                           inputs.E_s.yz.plan_f[n], inputs.E_s.yz.plan_b[n]);


          for (k = 1; k < K_tot; k++) {
            inputs.E_s.yz[k][j][i] =
                    PSTD.ca[n][k - 1] * inputs.E_s.yz[k][j][i] +
                    PSTD.cb[n][k - 1] * eh_vec[n][k][0] / ((double) PSTD.N_ez);
            // E_s.yz[k][j][i] = Enp1;
          }
        }
      // PSTD, E_s.yz
    }// if (solver_method == DerivativeMethod::FiniteDifference) (else
     // PseudoSpectral)
  }  // if(params.dimension==THREE || params.dimension==TE)

  if (inputs.params.dimension == THREE ||
      inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
    if (solver_method == SolverMethod::FiniteDifference) {
#pragma omp for
      // E_s.zx updates
      for (k = 0; k < K_tot; k++)
        for (j = 0; j < lv.J_tot_p1_bound; j++)
          for (i = 1; i < I_tot; i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (!inputs.params.is_multilayer) array_ind = i;
            else
              array_ind = (I_tot + 1) * k_loc + i;

            // use the average of material parameters between nodes
            if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
              rho = 0.;
              if (!inputs.materials[k][j][i]) {
                Ca = inputs.C.a.x[array_ind];
                Cb = inputs.C.b.x[array_ind];
                if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
                else
                  Cc = 0.;
              } else {
                Ca = inputs.Cmaterial.a.x[inputs.materials[k][j][i] - 1];
                Cb = inputs.Cmaterial.b.x[inputs.materials[k][j][i] - 1];
                Cc = inputs.Cmaterial.c.x[inputs.materials[k][j][i] - 1];
              }

              if (inputs.params.interp_mat_props) {
                if (!inputs.materials[k + 1][j][i]) {
                  Ca = Ca + inputs.C.a.x[array_ind];
                  Cb = Cb + inputs.C.b.x[array_ind];
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.x[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a.x[inputs.materials[k + 1][j][i] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.x[inputs.materials[k + 1][j][i] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.x[inputs.materials[k + 1][j][i] - 1];
                }

                Ca = Ca / 2.;
                Cb = Cb / 2.;
                Cc = Cc / 2.;
              }
            } else {
              Ca = inputs.C.a.x[array_ind];
              Cb = inputs.C.b.x[array_ind];
              if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
              else
                Cc = 0.;
              if (lv.is_conductive) rho = inputs.rho_cond.x[array_ind];
            }

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.x[array_ind];
              kappa_l = inputs.matched_layer.kappa.x[array_ind];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k + 1][j][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                  beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                  gamma_l += inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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
            Enp1 = Ca * inputs.E_s.zx[k][j][i] +
                   Cb * (inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i] -
                         inputs.H_s.yx[k][j][i - 1] -
                         inputs.H_s.yz[k][j][i - 1]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.zx[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dx *
                              ((1 + alpha_l) * lv.J_s.zx[k][j][i] +
                               beta_l * lv.J_nm1.zx[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dx * lv.J_c.zx[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.zx[k][j][i] +
                     beta_l * lv.J_nm1.zx[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.zx[k][j][i]);
              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
              lv.E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
              lv.J_nm1.zx[k][j][i] = lv.J_s.zx[k][j][i];
              lv.J_s.zx[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]);
            }

            inputs.E_s.zx[k][j][i] = Enp1;
          }
      // FDTD, E_s.zx
    } else {
#pragma omp for
      for (k = 0; k < K_tot; k++)
        for (j = 0; j < lv.J_tot_p1_bound; j++) {
          for (i = 1; i < I_tot; i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (!inputs.params.is_multilayer) array_ind = i;
            else
              array_ind = (I_tot + 1) * k_loc + i;

            // use the average of material parameters between nodes
            if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
              rho = 0.;
              if (!inputs.materials[k][j][i]) {
                Ca = inputs.C.a.x[array_ind];
                Cb = inputs.C.b.x[array_ind];
                if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
                else
                  Cc = 0.;
              } else {
                Ca = inputs.Cmaterial.a.x[inputs.materials[k][j][i] - 1];
                Cb = inputs.Cmaterial.b.x[inputs.materials[k][j][i] - 1];
                Cc = inputs.Cmaterial.c.x[inputs.materials[k][j][i] - 1];
              }

              if (inputs.params.interp_mat_props) {
                if (!inputs.materials[k + 1][j][i]) {
                  Ca = Ca + inputs.C.a.x[array_ind];
                  Cb = Cb + inputs.C.b.x[array_ind];
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.x[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a.x[inputs.materials[k + 1][j][i] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.x[inputs.materials[k + 1][j][i] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.x[inputs.materials[k + 1][j][i] - 1];
                }

                Ca = Ca / 2.;
                Cb = Cb / 2.;
                Cc = Cc / 2.;
              }
            } else {
              Ca = inputs.C.a.x[array_ind];
              Cb = inputs.C.b.x[array_ind];
              if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
              else
                Cc = 0.;
              if (lv.is_conductive) rho = inputs.rho_cond.x[array_ind];
            }

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.x[array_ind];
              kappa_l = inputs.matched_layer.kappa.x[array_ind];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k + 1][j][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                  beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                  gamma_l += inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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
            // Enp1 = Ca*E_s.zx[k][j][i]+Cb*(H_s.yx[k][j][i] + H_s.yz[k][j][i] -
            // H_s.yx[k][j][i-1] - H_s.yz[k][j][i-1]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.zx[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dx *
                              ((1 + alpha_l) * lv.J_s.zx[k][j][i] +
                               beta_l * lv.J_nm1.zx[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dx * lv.J_c.zx[k][j][i];
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.zx[k][j][i] +
                     beta_l * lv.J_nm1.zx[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.zx[k][j][i]);
              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
              lv.E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
              lv.J_nm1.zx[k][j][i] = lv.J_s.zx[k][j][i];
              lv.J_s.zx[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]);
            }

            eh_vec[n][i][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
            eh_vec[n][i][1] = 0.;
            PSTD.ca[n][i - 1] = Ca;
            PSTD.cb[n][i - 1] = Cb;
          }
          i = 0;
          eh_vec[n][i][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
          eh_vec[n][i][1] = 0.;

          first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ex, PSTD.N_ex,
                           inputs.E_s.zx.plan_f[n], inputs.E_s.zx.plan_b[n]);

          for (i = 1; i < I_tot; i++) {
            inputs.E_s.zx[k][j][i] =
                    PSTD.ca[n][i - 1] * inputs.E_s.zx[k][j][i] +
                    PSTD.cb[n][i - 1] * eh_vec[n][i][0] / ((double) PSTD.N_ex);
            // E_s.zx[k][j][i] = Enp1;
          }
        }
      // PSTD, E_s.zx
    }// if (solver_method == DerivativeMethod::FiniteDifference) (else
     // PseudoSpectral)
  }  //(params.dimension==THREE || params.dimension==TE)
  else {
#pragma omp for
    // E_s.zx updates
    for (k = 0; k <= K_tot; k++)
      for (j = 0; j < (J_tot + 1); j++)
        for (i = 1; i < I_tot; i++) {
          rho = 0.;
          k_loc = k;
          if (inputs.params.is_structure)
            if (k > inputs.params.pml.Dzl &&
                k < (inputs.params.pml.Dzl + lv.K)) {
              if ((k - inputs.structure[i][1]) <
                          (lv.K + inputs.params.pml.Dzl) &&
                  (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                k_loc = k - inputs.structure[i][1];
              else if ((k - inputs.structure[i][1]) >=
                       (lv.K + inputs.params.pml.Dzl))
                k_loc = inputs.params.pml.Dzl + lv.K - 1;
              else
                k_loc = inputs.params.pml.Dzl + 1;
            }
          if (!inputs.params.is_multilayer) array_ind = i;
          else
            array_ind = (I_tot + 1) * k_loc + i;

          // use the average of material parameters between nodes
          if (!inputs.materials[k][j][i]) {
            Ca = inputs.C.a.x[array_ind];
            Cb = inputs.C.b.x[array_ind];
            if (inputs.params.is_disp_ml) Cc = inputs.C.c.x[array_ind];
            else
              Cc = 0.;
            if (lv.is_conductive) rho = inputs.rho_cond.x[i];
          } else {
            rho = 0.;
            Ca = inputs.Cmaterial.a.x[inputs.materials[k][j][i] - 1];
            Cb = inputs.Cmaterial.b.x[inputs.materials[k][j][i] - 1];
            Cc = inputs.Cmaterial.c.x[inputs.materials[k][j][i] - 1];
          }

          alpha_l = 0.;
          beta_l = 0.;
          gamma_l = 0.;
          kappa_l = 1.;
          sigma_l = 0.;


          if (lv.is_disp || inputs.params.is_disp_ml) {
            sigma_l = inputs.matched_layer.sigma.x[array_ind];
            kappa_l = inputs.matched_layer.kappa.x[array_ind];
            alpha_l = inputs.matched_layer.alpha[k_loc];
            beta_l = inputs.matched_layer.beta[k_loc];
            gamma_l = inputs.matched_layer.gamma[k_loc];

            if (inputs.materials[k][j][i]) {
              alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
              beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
              gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];

            } else {
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
            }
          }

          Enp1 = Ca * inputs.E_s.zx[k][j][i] +
                 Cb * (inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i] -
                       inputs.H_s.yx[k][j][i - 1] - inputs.H_s.yz[k][j][i - 1]);
          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.zx[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dx *
                            ((1 + alpha_l) * lv.J_s.zx[k][j][i] +
                             beta_l * lv.J_nm1.zx[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dx * lv.J_c.zx[k][j][i];

          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
            Jnp1 = alpha_l * lv.J_s.zx[k][j][i] +
                   beta_l * lv.J_nm1.zx[k][j][i] +
                   kappa_l * gamma_l / (2. * inputs.params.dt) *
                           (Enp1 - lv.E_nm1.zx[k][j][i]);
            Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
            lv.E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
            lv.J_nm1.zx[k][j][i] = lv.J_s.zx[k][j][i];
            lv.J_s.zx[k][j][i] = Jnp1;
          }
          if (lv.is_conductive && rho) {
            lv.J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]);
          }

          inputs.E_s.zx[k][j][i] = Enp1;
        }
  }
  if (inputs.params.dimension == THREE ||
      inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
    if (solver_method == SolverMethod::FiniteDifference) {
      // FDTD, E_s.zy
#pragma omp for
      // E_s.zy updates
      for (k = 0; k < K_tot; k++)
        for (j = 1; j < J_tot; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (!inputs.params.is_multilayer) array_ind = j;
            else
              array_ind = (J_tot + 1) * k_loc + j;

            // use the average of material parameters between nodes
            if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
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
                if (!inputs.materials[k + 1][j][i]) {
                  Ca = Ca + inputs.C.a.y[array_ind];
                  Cb = Cb + inputs.C.b.y[array_ind];
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.y[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a.y[inputs.materials[k + 1][j][i] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.y[inputs.materials[k + 1][j][i] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.y[inputs.materials[k + 1][j][i] - 1];
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
                Cc = 0;
              if (lv.is_conductive) rho = inputs.rho_cond.y[array_ind];
            }

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.y[array_ind];
              kappa_l = inputs.matched_layer.kappa.y[array_ind];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k + 1][j][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                  beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                  gamma_l += inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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


            Enp1 = Ca * inputs.E_s.zy[k][j][i] +
                   Cb * (inputs.H_s.xy[k][j - 1][i] +
                         inputs.H_s.xz[k][j - 1][i] - inputs.H_s.xy[k][j][i] -
                         inputs.H_s.xz[k][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.zy[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dy *
                              ((1 + alpha_l) * lv.J_s.zy[k][j][i] +
                               beta_l * lv.J_nm1.zy[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dy * lv.J_c.zy[k][j][i];

            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.zy[k][j][i] +
                     beta_l * lv.J_nm1.zy[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.zy[k][j][i]);

              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
              lv.E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
              lv.J_nm1.zy[k][j][i] = lv.J_s.zy[k][j][i];
              lv.J_s.zy[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]);
            }
            inputs.E_s.zy[k][j][i] = Enp1;
          }
      // FDTD, E_s.zy
    } else {
#pragma omp for
      // E_s.zy updates
      for (k = 0; k < K_tot; k++)
        for (i = 0; i < (I_tot + 1); i++) {
          for (j = 1; j < J_tot; j++) {
            rho = 0.;
            k_loc = k;
            if (inputs.params.is_structure)
              if (k > inputs.params.pml.Dzl &&
                  k < (inputs.params.pml.Dzl + lv.K)) {
                if ((k - inputs.structure[i][1]) <
                            (lv.K + inputs.params.pml.Dzl) &&
                    (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                  k_loc = k - inputs.structure[i][1];
                else if ((k - inputs.structure[i][1]) >=
                         (lv.K + inputs.params.pml.Dzl))
                  k_loc = inputs.params.pml.Dzl + lv.K - 1;
                else
                  k_loc = inputs.params.pml.Dzl + 1;
              }
            if (!inputs.params.is_multilayer) array_ind = j;
            else
              array_ind = (J_tot + 1) * k_loc + j;

            // use the average of material parameters between nodes
            if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
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
                if (!inputs.materials[k + 1][j][i]) {
                  Ca = Ca + inputs.C.a.y[array_ind];
                  Cb = Cb + inputs.C.b.y[array_ind];
                  if (inputs.params.is_disp_ml)
                    Cc = Cc + inputs.C.c.y[array_ind];
                } else {
                  Ca = Ca +
                       inputs.Cmaterial.a.y[inputs.materials[k + 1][j][i] - 1];
                  Cb = Cb +
                       inputs.Cmaterial.b.y[inputs.materials[k + 1][j][i] - 1];
                  Cc = Cc +
                       inputs.Cmaterial.c.y[inputs.materials[k + 1][j][i] - 1];
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
                Cc = 0;
              if (lv.is_conductive) rho = inputs.rho_cond.y[array_ind];
            }

            alpha_l = 0.;
            beta_l = 0.;
            gamma_l = 0.;
            kappa_l = 1.;
            sigma_l = 0.;

            if (lv.is_disp || inputs.params.is_disp_ml) {
              sigma_l = inputs.matched_layer.sigma.y[array_ind];
              kappa_l = inputs.matched_layer.kappa.y[array_ind];
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
              if (inputs.materials[k][j][i] || inputs.materials[k + 1][j][i]) {
                if (inputs.materials[k][j][i]) {
                  alpha_l = inputs.alpha[inputs.materials[k][j][i] - 1];
                  beta_l = inputs.beta[inputs.materials[k][j][i] - 1];
                  gamma_l = inputs.gamma[inputs.materials[k][j][i] - 1];
                } else {
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                }

                if (inputs.materials[k + 1][j][i]) {
                  alpha_l += inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                  beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                  gamma_l += inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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


            // Enp1 = Ca*E_s.zy[k][j][i]+Cb*(H_s.xy[k][j-1][i] +
            // H_s.xz[k][j-1][i] - H_s.xy[k][j][i] - H_s.xz[k][j][i]);
            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
              Enp1 += Cc * lv.E_nm1.zy[k][j][i] -
                      1. / 2. * Cb * inputs.params.delta.dy *
                              ((1 + alpha_l) * lv.J_s.zy[k][j][i] +
                               beta_l * lv.J_nm1.zy[k][j][i]);
            if (lv.is_conductive && rho)
              Enp1 += Cb * inputs.params.delta.dy * lv.J_c.zy[k][j][i];

            if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
              Jnp1 = alpha_l * lv.J_s.zy[k][j][i] +
                     beta_l * lv.J_nm1.zy[k][j][i] +
                     kappa_l * gamma_l / (2. * inputs.params.dt) *
                             (Enp1 - lv.E_nm1.zy[k][j][i]);

              Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
              lv.E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
              lv.J_nm1.zy[k][j][i] = lv.J_s.zy[k][j][i];
              lv.J_s.zy[k][j][i] = Jnp1;
            }
            if (lv.is_conductive && rho) {
              lv.J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]);
            }

            eh_vec[n][j][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
            eh_vec[n][j][1] = 0.;
            PSTD.ca[n][j - 1] = Ca;
            PSTD.cb[n][j - 1] = Cb;
          }
          if (J_tot > 1) {
            j = 0;
            eh_vec[n][j][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
            eh_vec[n][j][1] = 0.;
            first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ey, PSTD.N_ey,
                             inputs.E_s.zy.plan_f[n], inputs.E_s.zy.plan_b[n]);
          }
          for (j = 1; j < J_tot; j++) {
            inputs.E_s.zy[k][j][i] =
                    PSTD.ca[n][j - 1] * inputs.E_s.zy[k][j][i] -
                    PSTD.cb[n][j - 1] * eh_vec[n][j][0] / ((double) PSTD.N_ey);
            // E_s.zy[k][j][i] = Enp1;
          }
        }
      // PSTD, E_s.zy
    }// if (solver_method == DerivativeMethod::FiniteDifference) (else
     // PseudoSpectral)
  }  //(params.dimension==THREE || params.dimension==TE)
  else {
#pragma omp for
    for (k = 0; k <= K_tot; k++)
      for (j = 1; j < J_tot; j++)
        for (i = 0; i < (I_tot + 1); i++) {
          rho = 0.;
          k_loc = k;
          if (inputs.params.is_structure)
            if (k > inputs.params.pml.Dzl &&
                k < (inputs.params.pml.Dzl + lv.K)) {
              if ((k - inputs.structure[i][1]) <
                          (lv.K + inputs.params.pml.Dzl) &&
                  (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                k_loc = k - inputs.structure[i][1];
              else if ((k - inputs.structure[i][1]) >=
                       (lv.K + inputs.params.pml.Dzl))
                k_loc = inputs.params.pml.Dzl + lv.K - 1;
              else
                k_loc = inputs.params.pml.Dzl + 1;
            }
          if (!inputs.params.is_multilayer) array_ind = j;
          else
            array_ind = (J_tot + 1) * k_loc + j;

          // use the average of material parameters between nodes
          if (!inputs.materials[k][j][i]) {
            Ca = inputs.C.a.y[array_ind];
            Cb = inputs.C.b.y[array_ind];
            if (inputs.params.is_disp_ml) Cc = inputs.C.c.y[array_ind];
            else
              Cc = 0.;
            if (lv.is_conductive) rho = inputs.rho_cond.y[array_ind];
          } else {
            rho = 0.;
            Ca = inputs.Cmaterial.a.y[inputs.materials[k][j][i] - 1];
            Cb = inputs.Cmaterial.b.y[inputs.materials[k][j][i] - 1];
            Cc = inputs.Cmaterial.c.y[inputs.materials[k][j][i] - 1];
          }

          alpha_l = 0.;
          beta_l = 0.;
          gamma_l = 0.;
          kappa_l = 1.;
          sigma_l = 0.;

          if (lv.is_disp || inputs.params.is_disp_ml) {
            kappa_l = inputs.matched_layer.kappa.y[array_ind];
            sigma_l = inputs.matched_layer.sigma.y[array_ind];
            alpha_l = inputs.matched_layer.alpha[k_loc];
            beta_l = inputs.matched_layer.beta[k_loc];
            gamma_l = inputs.matched_layer.gamma[k_loc];

            if (!inputs.materials[k][j][i]) {
              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
            } else {
              alpha_l = inputs.matched_layer.alpha[k_loc];
              beta_l = inputs.matched_layer.beta[k_loc];
              gamma_l = inputs.matched_layer.gamma[k_loc];
            }
          }


          Enp1 = Ca * inputs.E_s.zy[k][j][i] +
                 Cb * (inputs.H_s.xy[k][j - 1][i] + inputs.H_s.xz[k][j - 1][i] -
                       inputs.H_s.xy[k][j][i] - inputs.H_s.xz[k][j][i]);
          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.zy[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dy *
                            ((1 + alpha_l) * lv.J_s.zy[k][j][i] +
                             beta_l * lv.J_nm1.zy[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dy * lv.J_c.zy[k][j][i];

          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
            Jnp1 = alpha_l * lv.J_s.zy[k][j][i] +
                   beta_l * lv.J_nm1.zy[k][j][i] +
                   kappa_l * gamma_l / (2. * inputs.params.dt) *
                           (Enp1 - lv.E_nm1.zy[k][j][i]);

            Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
            lv.E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
            lv.J_nm1.zy[k][j][i] = lv.J_s.zy[k][j][i];
            lv.J_s.zy[k][j][i] = Jnp1;
          }
          if (lv.is_conductive && rho) {
            lv.J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]);
          }

          inputs.E_s.zy[k][j][i] = Enp1;
        }
  }// end of parallel section
}

void SimulationManager::update_Exy(LoopVariables &lv) {
  // Fetch simulation dimensions for easy asscess
  int I_tot = n_Yee_cells().I_tot(), J_tot = n_Yee_cells().J_tot(),
      K_tot = n_Yee_cells().K_tot();
  // Get the thread number
  int n = omp_get_thread_num();
#pragma omp for
  for (int k = 0; k < (K_tot + 1); k++)
    for (int i = 0; i < I_tot; i++) {
      for (int j = 1; j < J_tot; j++) {
        double rho = 0.;
        int k_loc = k, array_ind;
        double Ca, Cb, Cc;

        if (inputs.params.is_structure) {
          if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + lv.K)) {
            if ((k - inputs.structure[i][1]) < (lv.K + inputs.params.pml.Dzl) &&
                (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
              k_loc = k - inputs.structure[i][1];
            else if ((k - inputs.structure[i][1]) >=
                     (lv.K + inputs.params.pml.Dzl))
              k_loc = inputs.params.pml.Dzl + lv.K - 1;
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

        if (lv.is_disp || inputs.params.is_disp_ml) {
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
        // DIFFERENCES!
        double Enp1, Jnp1;
        if (solver_method == SolverMethod::FiniteDifference) {
          Enp1 = Ca * inputs.E_s.xy[k][j][i] +
                 Cb * (inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i] -
                       inputs.H_s.zy[k][j - 1][i] - inputs.H_s.zx[k][j - 1][i]);
          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.xy[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dy *
                            ((1 + alpha_l) * lv.J_s.xy[k][j][i] +
                             beta_l * lv.J_nm1.xy[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dy * lv.J_c.xy[k][j][i];
          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
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
          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l)
            Enp1 += Cc * lv.E_nm1.xy[k][j][i] -
                    1. / 2. * Cb * inputs.params.delta.dy *
                            ((1 + alpha_l) * lv.J_s.xy[k][j][i] +
                             beta_l * lv.J_nm1.xy[k][j][i]);
          if (lv.is_conductive && rho)
            Enp1 += Cb * inputs.params.delta.dy * lv.J_c.xy[k][j][i];
          if ((lv.is_disp || inputs.params.is_disp_ml) && gamma_l) {
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
