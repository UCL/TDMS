/**
 * @file execute_simulation.cpp
 * @brief Code for performing the time-stepping algorithm for the physical
 * process TDMS models.
 */
#include "simulation_manager/simulation_manager.h"

#include <omp.h>
#include <spdlog/spdlog.h>

#include "numerical_derivative.h"
#include "simulation_manager/loop_variables.h"

using namespace std;
using namespace tdms_math_constants;
using namespace tdms_phys_constants;
using tdms_flags::SolverMethod;

void SimulationManager::execute() {
  // log the number of OMP threads being used
  spdlog::info("Using {} OMP threads", omp_get_max_threads());

  // get the number of Yee cells in each axial direction
  IJKDimensions IJK_tot = n_Yee_cells();
  int I_tot = IJK_tot.i, J_tot = IJK_tot.j, K_tot = IJK_tot.k;

  // DECLARE VARIABLES SCOPED TO THIS FUNCTION ONLY
  double rho;
  double alpha_l, beta_l, gamma_l;//< alpha, beta, gamma parameters of the layer
                                  // the local thread is examining
  double kappa_l, sigma_l;//< kappa, sigma parameters of the layer the local
                          // thread is examining
  double time_of_last_log_write;//< (Real) time since the last iteration log was
                                // written to the screen

  double Ca, Cb, Cc;// used by interpolation scheme
  // the C and D vars for free space and pml
  double Enp1, Jnp1;

  int i, j, k;//< Loop variables
  int n;      //< The thread number of the local OMP thread
  int k_loc;  //< Local thread copy of the variable k

  int dft_counter = 0;//< Number of DFTs we have performed since last checking
                      // for phasor convergence

  complex<double> Idxt, Idyt, kprop;

  // variables used in the main loop that require linking/setup from the input
  // and output objects
  LoopVariables loop_variables(inputs, outputs.get_E_dimensions());

  /*The times of the E and H fields at the point where update equations are
    applied. time_H is actually the time of the H field when the E field
    consistency update is applied and vice versa. time_E > time_H below since
    after the E field consistency update the E field will have advanced one time
    step.

    The interpretation of time is slightly complicated in the following. In what
    follows I write (tind*dt,(tind+1/2)*dt) to mean the time at which we
    currently know the electric (tind*dt) and magnetic ( (tind+1/2)*dt ) fields.

    Times before                Operation         Times after
    (tind*dt,(tind+1/2)*dt)     Extract phasors   (tind*dt,(tind+1/2)*dt)
    (tind*dt,(tind+1/2)*dt)     E field update    ( (tind+1)*dt,(tind+1/2)*dt)
    ((tind+1)*dt,(tind+1/2)*dt) H field update    ( (tind+1)*dt,(tind+3/2)*dt)
    ((tind+1)*dt,(tind+3/2)*dt) Normalisation extraction

    We note that the extractPhasorENorm uses (tind+1)*dt and extractPhasorHNorm
    uses (tind+1/2)*dt to perform the update equation in the DFT. This seems
    incorrect at first but we note that they take the terms fte and fth as
    inputs respectively. When one notes that fte is calculated using time_E and
    fth using time_H we see that this indexing is correct, ie, time_E =
    (tind+1)*dt and time_H = (tind+1/2)*dt.
  */
  double time_E, time_H;

  // fetch the current time for logging purposes
  time_of_last_log_write = (double) time(NULL);
  spdlog::info("Starting main loop");

  if (TIME_MAIN_LOOP) { timers.start_timer(TimersTrackingLoop::MAIN); }

  for (unsigned int tind = inputs.params.start_tind; tind < inputs.params.Nt;
       tind++) {
    time_E = ((double) (tind + 1)) * inputs.params.dt;
    time_H = time_E - inputs.params.dt / 2.;

    // Check for phasor convergence, break if achieved
    if (check_phasor_convergence(dft_counter,
                                 loop_variables.E_at_previous_iteration)) {
      break;
    }

    // Extract the volume, surface, and vertex phasors to the output
    extract_phasors(dft_counter, tind);

    // Extract the fields at the sample locations
    if (outputs.fieldsample.all_vectors_are_non_empty()) {
      outputs.fieldsample.extract(inputs.E_s, inputs.params.pml,
                                  inputs.params.Nt);
    }

    // Compute the detector function
    compute_detector_functions(tind, loop_variables);

    // Update equations for the E field

    /*There are two options for determining the update coefficients for the FDTD
      cell:

      1) If cell (i,j,k) is either free space or PML:

      materials[k][j][i] will be set to 0. In this case the update parameter
      used will be given by C.a.y[j], C.b.y[j] etc depending on which update
      equation is being implemented.

      2) if cell (i,j,k) is composed of a scattering type material then
      materials[k][j][i] will be non-zero and will be an index into
      Cmaterial.a.y and Cmaterial.b.y etc depending on which update equation is
      being implemented.

    */

    int array_ind = 0;
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
    (void) n;// n is unused in FD derivatives – this silences the compiler
             // warning

#pragma omp parallel default(shared) private(                                  \
        i, j, k, n, rho, k_loc, array_ind, Ca, Cb, Cc, alpha_l, beta_l,        \
        gamma_l, kappa_l, sigma_l, Enp1, Jnp1)//,ca_vec,cb_vec,eh_vec)
    {
      n = omp_get_thread_num();
      Enp1 = 0.0;
      array_ind = 0;

      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        update_Exy(loop_variables);
        update_Exz(loop_variables);
        // E_s.yx updates
        if (solver_method == SolverMethod::FiniteDifference) {
          // FDTD, E_s.yx
#pragma omp for
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound; j++)
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure) {
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
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
                      Ca = Ca + inputs.Cmaterial.a.x[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cb = Cb + inputs.Cmaterial.b.x[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cc = Cc + inputs.Cmaterial.c.x[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
                                              1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot,
                                                                    j + 1)][i] -
                                            1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
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


                Enp1 = Ca * inputs.E_s.yx(i, j, k) +
                       Cb * (inputs.H_s.zx(i - 1, j, k) +
                             inputs.H_s.zy(i - 1, j, k) -
                             inputs.H_s.zx(i, j, k) - inputs.H_s.zy(i, j, k));
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yx(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.yx(i, j, k) +
                                   beta_l * loop_variables.J_nm1.yx(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx *
                          loop_variables.J_c.yx(i, j, k);
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yx(i, j, k) +
                         beta_l * loop_variables.J_nm1.yx(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yx(i, j, k));
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx(i, j, k);
                  loop_variables.E_nm1.yx(i, j, k) = inputs.E_s.yx(i, j, k);
                  loop_variables.J_nm1.yx(i, j, k) =
                          loop_variables.J_s.yx(i, j, k);
                  loop_variables.J_s.yx(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yx(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.yx(i, j, k));
                }

                inputs.E_s.yx(i, j, k) = Enp1;
              }
          // FDTD, E_s.yx
        } else {
#pragma omp for
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound; j++) {
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure) {
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
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
                      Ca = Ca + inputs.Cmaterial.a.x[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cb = Cb + inputs.Cmaterial.b.x[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cc = Cc + inputs.Cmaterial.c.x[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
                                              1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot,
                                                                    j + 1)][i] -
                                            1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
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


                // Enp1 = Ca*E_s.yx(i,j,k)+Cb*(H_s.zx[k][j][i-1] +
                // H_s.zy[k][j][i-1] - H_s.zx(i,j,k) - H_s.zy(i,j,k));
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yx(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.yx(i, j, k) +
                                   beta_l * loop_variables.J_nm1.yx(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx *
                          loop_variables.J_c.yx(i, j, k);
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yx(i, j, k) +
                         beta_l * loop_variables.J_nm1.yx(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yx(i, j, k));
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx(i, j, k);
                  loop_variables.E_nm1.yx(i, j, k) = inputs.E_s.yx(i, j, k);
                  loop_variables.J_nm1.yx(i, j, k) =
                          loop_variables.J_s.yx(i, j, k);
                  loop_variables.J_s.yx(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yx(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.yx(i, j, k));
                }

                eh_vec[n][i][0] =
                        inputs.H_s.zx(i, j, k) + inputs.H_s.zy(i, j, k);
                eh_vec[n][i][1] = 0.;
                PSTD.ca(n, i - 1) = Ca;
                PSTD.cb(n, i - 1) = Cb;
              }
              i = 0;
              eh_vec[n][i][0] = inputs.H_s.zx(i, j, k) + inputs.H_s.zy(i, j, k);
              eh_vec[n][i][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ex, PSTD.N_ex,
                               inputs.E_s.yx.plan_f[n],
                               inputs.E_s.yx.plan_b[n]);

              for (i = 1; i < I_tot; i++) {
                inputs.E_s.yx(i, j, k) =
                        PSTD.ca(n, i - 1) * inputs.E_s.yx(i, j, k) -
                        PSTD.cb(n, i - 1) * eh_vec[n][i][0] /
                                ((double) PSTD.N_ex);
                // E_s.yx(i,j,k) = Enp1;
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
            for (j = 0; j < loop_variables.J_loop_upper_bound; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
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
                      if (inputs.params.is_disp_ml)
                        Cc = Cc + inputs.C.c.z[k_loc];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.z[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cb = Cb + inputs.Cmaterial.b.z[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cc = Cc + inputs.Cmaterial.c.z[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
                                              1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot,
                                                                    j + 1)][i] -
                                            1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
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
                Enp1 = Ca * inputs.E_s.yz(i, j, k) +
                       Cb * (inputs.H_s.xy(i, j, k) + inputs.H_s.xz(i, j, k) -
                             inputs.H_s.xy(i, j, k - 1) -
                             inputs.H_s.xz(i, j, k - 1));
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yz(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.yz(i, j, k) +
                                   beta_l * loop_variables.J_nm1.yz(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dz *
                          loop_variables.J_c.yz(i, j, k);

                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yz(i, j, k) +
                         beta_l * loop_variables.J_nm1.yz(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yz(i, j, k));
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz(i, j, k);
                  loop_variables.E_nm1.yz(i, j, k) = inputs.E_s.yz(i, j, k);
                  loop_variables.J_nm1.yz(i, j, k) =
                          loop_variables.J_s.yz(i, j, k);
                  loop_variables.J_s.yz(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yz(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.yz(i, j, k));
                }

                inputs.E_s.yz(i, j, k) = Enp1;
              }
          // FDTD, E_s.yz
        } else {
#pragma omp for
          for (j = 0; j < loop_variables.J_loop_upper_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (k = 1; k < K_tot; k++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
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
                      if (inputs.params.is_disp_ml)
                        Cc = Cc + inputs.C.c.z[k_loc];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.z[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cb = Cb + inputs.Cmaterial.b.z[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
                                                     1];
                      Cc = Cc + inputs.Cmaterial.c.z[inputs.materials[k][min(
                                                             J_tot, j + 1)][i] -
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
                                              1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot,
                                                                    j + 1)][i] -
                                            1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(
                                                      J_tot, j + 1)][i] -
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
                // Enp1 = Ca*E_s.yz(i,j,k)+Cb*(H_s.xy(i,j,k) +
                // H_s.xz(i, j, k) - H_s.xy[k-1][j][i] - H_s.xz[k-1][j][i]);
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yz(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.yz(i, j, k) +
                                   beta_l * loop_variables.J_nm1.yz(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dz *
                          loop_variables.J_c.yz(i, j, k);

                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yz(i, j, k) +
                         beta_l * loop_variables.J_nm1.yz(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yz(i, j, k));
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz(i, j, k);
                  loop_variables.E_nm1.yz(i, j, k) = inputs.E_s.yz(i, j, k);
                  loop_variables.J_nm1.yz(i, j, k) =
                          loop_variables.J_s.yz(i, j, k);
                  loop_variables.J_s.yz(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yz(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.yz(i, j, k));
                }

                eh_vec[n][k][0] =
                        inputs.H_s.xy(i, j, k) + inputs.H_s.xz(i, j, k);
                eh_vec[n][k][1] = 0.;
                PSTD.ca(n, k - 1) = Ca;
                PSTD.cb(n, k - 1) = Cb;
              }
              k = 0;
              eh_vec[n][k][0] = inputs.H_s.xy(i, j, k) + inputs.H_s.xz(i, j, k);
              eh_vec[n][k][1] = 0.;
              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ez, PSTD.N_ez,
                               inputs.E_s.yz.plan_f[n],
                               inputs.E_s.yz.plan_b[n]);


              for (k = 1; k < K_tot; k++) {
                inputs.E_s.yz(i, j, k) =
                        PSTD.ca(n, k - 1) * inputs.E_s.yz(i, j, k) +
                        PSTD.cb(n, k - 1) * eh_vec[n][k][0] /
                                ((double) PSTD.N_ez);
                // E_s.yz(i,j,k) = Enp1;
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
            for (j = 0; j < loop_variables.J_loop_upper_bound_plus_1; j++)
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;

                // use the average of material parameters between nodes
                if (inputs.materials[k][j][i] ||
                    inputs.materials[k + 1][j][i]) {
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
                      Ca = Ca + inputs.Cmaterial.a
                                        .x[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b
                                        .x[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c
                                        .x[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.x[array_ind];
                  kappa_l = inputs.matched_layer.kappa.x[array_ind];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] ||
                      inputs.materials[k + 1][j][i]) {
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
                      alpha_l +=
                              inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                      gamma_l +=
                              inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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
                Enp1 = Ca * inputs.E_s.zx(i, j, k) +
                       Cb * (inputs.H_s.yx(i, j, k) + inputs.H_s.yz(i, j, k) -
                             inputs.H_s.yx(i - 1, j, k) -
                             inputs.H_s.yz(i - 1, j, k));
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zx(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.zx(i, j, k) +
                                   beta_l * loop_variables.J_nm1.zx(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx *
                          loop_variables.J_c.zx(i, j, k);
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zx(i, j, k) +
                         beta_l * loop_variables.J_nm1.zx(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zx(i, j, k));
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx(i, j, k);
                  loop_variables.E_nm1.zx(i, j, k) = inputs.E_s.zx(i, j, k);
                  loop_variables.J_nm1.zx(i, j, k) =
                          loop_variables.J_s.zx(i, j, k);
                  loop_variables.J_s.zx(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zx(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.zx(i, j, k));
                }

                inputs.E_s.zx(i, j, k) = Enp1;
              }
          // FDTD, E_s.zx
        } else {
#pragma omp for
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound_plus_1; j++) {
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;

                // use the average of material parameters between nodes
                if (inputs.materials[k][j][i] ||
                    inputs.materials[k + 1][j][i]) {
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
                      Ca = Ca + inputs.Cmaterial.a
                                        .x[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b
                                        .x[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c
                                        .x[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.x[array_ind];
                  kappa_l = inputs.matched_layer.kappa.x[array_ind];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] ||
                      inputs.materials[k + 1][j][i]) {
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
                      alpha_l +=
                              inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                      gamma_l +=
                              inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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
                // Enp1 = Ca*E_s.zx(i,j,k)+Cb*(H_s.yx(i, j, k) +
                // H_s.yz(i,j,k) - H_s.yx[k][j][i-1] - H_s.yz[k][j][i-1]);
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zx(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.zx(i, j, k) +
                                   beta_l * loop_variables.J_nm1.zx(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx *
                          loop_variables.J_c.zx(i, j, k);
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zx(i, j, k) +
                         beta_l * loop_variables.J_nm1.zx(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zx(i, j, k));
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx(i, j, k);
                  loop_variables.E_nm1.zx(i, j, k) = inputs.E_s.zx(i, j, k);
                  loop_variables.J_nm1.zx(i, j, k) =
                          loop_variables.J_s.zx(i, j, k);
                  loop_variables.J_s.zx(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zx(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.zx(i, j, k));
                }

                eh_vec[n][i][0] =
                        inputs.H_s.yx(i, j, k) + inputs.H_s.yz(i, j, k);
                eh_vec[n][i][1] = 0.;
                PSTD.ca(n, i - 1) = Ca;
                PSTD.cb(n, i - 1) = Cb;
              }
              i = 0;
              eh_vec[n][i][0] = inputs.H_s.yx(i, j, k) + inputs.H_s.yz(i, j, k);
              eh_vec[n][i][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ex, PSTD.N_ex,
                               inputs.E_s.zx.plan_f[n],
                               inputs.E_s.zx.plan_b[n]);

              for (i = 1; i < I_tot; i++) {
                inputs.E_s.zx(i, j, k) =
                        PSTD.ca(n, i - 1) * inputs.E_s.zx(i, j, k) +
                        PSTD.cb(n, i - 1) * eh_vec[n][i][0] /
                                ((double) PSTD.N_ex);
                // E_s.zx(i,j,k) = Enp1;
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
                    k < (inputs.params.pml.Dzl +
                         loop_variables.n_non_pml_cells_in_K)) {
                  if ((k - inputs.structure[i][1]) <
                              (loop_variables.n_non_pml_cells_in_K +
                               inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.n_non_pml_cells_in_K +
                            inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl +
                            loop_variables.n_non_pml_cells_in_K - 1;
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
                if (loop_variables.is_conductive) rho = inputs.rho_cond.x[i];
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


              if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
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

              Enp1 = Ca * inputs.E_s.zx(i, j, k) +
                     Cb * (inputs.H_s.yx(i, j, k) + inputs.H_s.yz(i, j, k) -
                           inputs.H_s.yx(i - 1, j, k) -
                           inputs.H_s.yz(i - 1, j, k));
              if ((loop_variables.is_dispersive || inputs.params.is_disp_ml) &&
                  gamma_l)
                Enp1 += Cc * loop_variables.E_nm1.zx(i, j, k) -
                        1. / 2. * Cb * inputs.params.delta.dx *
                                ((1 + alpha_l) *
                                         loop_variables.J_s.zx(i, j, k) +
                                 beta_l * loop_variables.J_nm1.zx(i, j, k));
              if (loop_variables.is_conductive && rho)
                Enp1 += Cb * inputs.params.delta.dx *
                        loop_variables.J_c.zx(i, j, k);

              if ((loop_variables.is_dispersive || inputs.params.is_disp_ml) &&
                  gamma_l) {
                Jnp1 = alpha_l * loop_variables.J_s.zx(i, j, k) +
                       beta_l * loop_variables.J_nm1.zx(i, j, k) +
                       kappa_l * gamma_l / (2. * inputs.params.dt) *
                               (Enp1 - loop_variables.E_nm1.zx(i, j, k));
                Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx(i, j, k);
                loop_variables.E_nm1.zx(i, j, k) = inputs.E_s.zx(i, j, k);
                loop_variables.J_nm1.zx(i, j, k) =
                        loop_variables.J_s.zx(i, j, k);
                loop_variables.J_s.zx(i, j, k) = Jnp1;
              }
              if (loop_variables.is_conductive && rho) {
                loop_variables.J_c.zx(i, j, k) -=
                        rho * (Enp1 + inputs.E_s.zx(i, j, k));
              }

              inputs.E_s.zx(i, j, k) = Enp1;
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
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;

                // use the average of material parameters between nodes
                if (inputs.materials[k][j][i] ||
                    inputs.materials[k + 1][j][i]) {
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
                      Ca = Ca + inputs.Cmaterial.a
                                        .y[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b
                                        .y[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c
                                        .y[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.y[array_ind];
                  kappa_l = inputs.matched_layer.kappa.y[array_ind];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] ||
                      inputs.materials[k + 1][j][i]) {
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
                      alpha_l +=
                              inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                      gamma_l +=
                              inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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


                Enp1 = Ca * inputs.E_s.zy(i, j, k) +
                       Cb * (inputs.H_s.xy(i, j - 1, k) +
                             inputs.H_s.xz(i, j - 1, k) -
                             inputs.H_s.xy(i, j, k) - inputs.H_s.xz(i, j, k));
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zy(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.zy(i, j, k) +
                                   beta_l * loop_variables.J_nm1.zy(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dy *
                          loop_variables.J_c.zy(i, j, k);

                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zy(i, j, k) +
                         beta_l * loop_variables.J_nm1.zy(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zy(i, j, k));

                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy(i, j, k);
                  loop_variables.E_nm1.zy(i, j, k) = inputs.E_s.zy(i, j, k);
                  loop_variables.J_nm1.zy(i, j, k) =
                          loop_variables.J_s.zy(i, j, k);
                  loop_variables.J_s.zy(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zy(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.zy(i, j, k));
                }
                inputs.E_s.zy(i, j, k) = Enp1;
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
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;

                // use the average of material parameters between nodes
                if (inputs.materials[k][j][i] ||
                    inputs.materials[k + 1][j][i]) {
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
                      Ca = Ca + inputs.Cmaterial.a
                                        .y[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b
                                        .y[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c
                                        .y[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive)
                    rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.y[array_ind];
                  kappa_l = inputs.matched_layer.kappa.y[array_ind];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] ||
                      inputs.materials[k + 1][j][i]) {
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
                      alpha_l +=
                              inputs.alpha[inputs.materials[k + 1][j][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k + 1][j][i] - 1];
                      gamma_l +=
                              inputs.gamma[inputs.materials[k + 1][j][i] - 1];
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


                // Enp1 = Ca*E_s.zy(i,j,k)+Cb*(H_s.xy[k][j-1][i] +
                // H_s.xz[k][j-1][i] - H_s.xy(i,j,k) - H_s.xz(i, j, k));
                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zy(i, j, k) -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) *
                                           loop_variables.J_s.zy(i, j, k) +
                                   beta_l * loop_variables.J_nm1.zy(i, j, k));
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dy *
                          loop_variables.J_c.zy(i, j, k);

                if ((loop_variables.is_dispersive ||
                     inputs.params.is_disp_ml) &&
                    gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zy(i, j, k) +
                         beta_l * loop_variables.J_nm1.zy(i, j, k) +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zy(i, j, k));

                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy(i, j, k);
                  loop_variables.E_nm1.zy(i, j, k) = inputs.E_s.zy(i, j, k);
                  loop_variables.J_nm1.zy(i, j, k) =
                          loop_variables.J_s.zy(i, j, k);
                  loop_variables.J_s.zy(i, j, k) = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zy(i, j, k) -=
                          rho * (Enp1 + inputs.E_s.zy(i, j, k));
                }

                eh_vec[n][j][0] =
                        inputs.H_s.xy(i, j, k) + inputs.H_s.xz(i, j, k);
                eh_vec[n][j][1] = 0.;
                PSTD.ca(n, j - 1) = Ca;
                PSTD.cb(n, j - 1) = Cb;
              }
              if (J_tot > 1) {
                j = 0;
                eh_vec[n][j][0] =
                        inputs.H_s.xy(i, j, k) + inputs.H_s.xz(i, j, k);
                eh_vec[n][j][1] = 0.;
                first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ey, PSTD.N_ey,
                                 inputs.E_s.zy.plan_f[n],
                                 inputs.E_s.zy.plan_b[n]);
              }
              for (j = 1; j < J_tot; j++) {
                inputs.E_s.zy(i, j, k) =
                        PSTD.ca(n, j - 1) * inputs.E_s.zy(i, j, k) -
                        PSTD.cb(n, j - 1) * eh_vec[n][j][0] /
                                ((double) PSTD.N_ey);
                // E_s.zy(i,j,k) = Enp1;
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
                    k < (inputs.params.pml.Dzl +
                         loop_variables.n_non_pml_cells_in_K)) {
                  if ((k - inputs.structure[i][1]) <
                              (loop_variables.n_non_pml_cells_in_K +
                               inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.n_non_pml_cells_in_K +
                            inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl +
                            loop_variables.n_non_pml_cells_in_K - 1;
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
                if (loop_variables.is_conductive)
                  rho = inputs.rho_cond.y[array_ind];
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

              if (loop_variables.is_dispersive || inputs.params.is_disp_ml) {
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


              Enp1 = Ca * inputs.E_s.zy(i, j, k) +
                     Cb * (inputs.H_s.xy(i, j - 1, k) +
                           inputs.H_s.xz(i, j - 1, k) - inputs.H_s.xy(i, j, k) -
                           inputs.H_s.xz(i, j, k));
              if ((loop_variables.is_dispersive || inputs.params.is_disp_ml) &&
                  gamma_l)
                Enp1 += Cc * loop_variables.E_nm1.zy(i, j, k) -
                        1. / 2. * Cb * inputs.params.delta.dy *
                                ((1 + alpha_l) *
                                         loop_variables.J_s.zy(i, j, k) +
                                 beta_l * loop_variables.J_nm1.zy(i, j, k));
              if (loop_variables.is_conductive && rho)
                Enp1 += Cb * inputs.params.delta.dy *
                        loop_variables.J_c.zy(i, j, k);

              if ((loop_variables.is_dispersive || inputs.params.is_disp_ml) &&
                  gamma_l) {
                Jnp1 = alpha_l * loop_variables.J_s.zy(i, j, k) +
                       beta_l * loop_variables.J_nm1.zy(i, j, k) +
                       kappa_l * gamma_l / (2. * inputs.params.dt) *
                               (Enp1 - loop_variables.E_nm1.zy(i, j, k));

                Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy(i, j, k);
                loop_variables.E_nm1.zy(i, j, k) = inputs.E_s.zy(i, j, k);
                loop_variables.J_nm1.zy(i, j, k) =
                        loop_variables.J_s.zy(i, j, k);
                loop_variables.J_s.zy(i, j, k) = Jnp1;
              }
              if (loop_variables.is_conductive && rho) {
                loop_variables.J_c.zy(i, j, k) -=
                        rho * (Enp1 + inputs.E_s.zy(i, j, k));
              }

              inputs.E_s.zy(i, j, k) = Enp1;
            }
      }
    }// end of parallel section
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
    /********************/

    /* Update source terms for self consistency across scattered/total
     * interface. E_s updates use time_H simulation time. */

    if (inputs.params.source_mode == SourceMode::steadystate) {
      // Steady-state source term updates
      E_source_update_all_steadystate(time_H, loop_variables.is_conductive,
                                      loop_variables.J_c, loop_variables.J_s);

      // Common phase term in update equations
      complex<double> common_phase =
              exp(-IMAGINARY_UNIT *
                  fmod(inputs.params.omega_an * time_H, 2. * DCPI));
      // Common amplitude factor in update equations
      double common_amplitude = linear_ramp(time_H);
      // Update output H-field
      outputs.H.ft = real(common_amplitude * common_phase);
    } else if (inputs.params.source_mode == SourceMode::pulsed) {
      // Pulsed source term updates
      update_source_terms_pulsed(time_H, loop_variables.is_conductive,
                                 loop_variables.J_c, loop_variables.J_s);

      // Common amplitude factor in update equations
      double common_amplitude =
              exp(-1.0 * DCPI *
                  pow((time_H - inputs.params.to_l +
                       inputs.params.delta.dz / LIGHT_V / 2.) /
                              inputs.params.hwhm,
                      2));
      // Common phase term in update equations
      complex<double> common_phase =
              -1.0 * IMAGINARY_UNIT *
              exp(-IMAGINARY_UNIT *
                  fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                       2. * DCPI));
      // Update output H-field
      outputs.H.ft = real(common_phase) * common_amplitude;
    }

    // end of source terms
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    /********************/
    // begin parallel
#pragma omp parallel default(shared) private(                                  \
        i, j, k, n, k_loc, array_ind)//,ca_vec,cb_vec,eh_vec)
    {
      n = omp_get_thread_num();

      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, H_s.xz
#pragma omp for
          // H_s.xz updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }

                if (!inputs.materials[k][j][i])
                  inputs.H_s.xz(i, j, k) =
                          inputs.D.a.z[k_loc] * inputs.H_s.xz(i, j, k) +
                          inputs.D.b.z[k_loc] * (inputs.E_s.yx(i, j, k + 1) +
                                                 inputs.E_s.yz(i, j, k + 1) -
                                                 inputs.E_s.yx(i, j, k) -
                                                 inputs.E_s.yz(i, j, k));
                else
                  inputs.H_s.xz(i, j, k) =
                          inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.xz(i, j, k) +
                          inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.yx(i, j, k + 1) +
                                   inputs.E_s.yz(i, j, k + 1) -
                                   inputs.E_s.yx(i, j, k) -
                                   inputs.E_s.yz(i, j, k));
              }
          // FDTD, H_s.xz
        } else {
#pragma omp for
          // H_s.xz updates
          for (j = 0; j < loop_variables.J_loop_upper_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }

                if (!inputs.materials[k][j][i]) {
                  PSTD.ca(n, k) = inputs.D.a.z[k_loc];
                  PSTD.cb(n, k) = inputs.D.b.z[k_loc];
                } else {
                  PSTD.ca(n, k) =
                          inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1];
                  PSTD.cb(n, k) =
                          inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1];
                  // H_s.xz(i, j, k) =
                  // Dmaterial.Da.z[materials[k][j][i]-1]*H_s.xz(i, j,
                  // k)+Dmaterial.Db.z[materials[k][j][i]-1]*(E_s.yx[k+1][j][i]
                  // + E_s.yz[k+1][j][i] - E_s.yx(i,j,k) - E_s.yz(i,j,k));
                }

                eh_vec[n][k][0] =
                        inputs.E_s.yx(i, j, k) + inputs.E_s.yz(i, j, k);
                eh_vec[n][k][1] = 0.;
              }
              k = K_tot;
              eh_vec[n][k][0] = inputs.E_s.yx(i, j, k) + inputs.E_s.yz(i, j, k);
              eh_vec[n][k][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hz, PSTD.N_hz,
                               inputs.H_s.xz.plan_f[n],
                               inputs.H_s.xz.plan_b[n]);

              for (k = 0; k < K_tot; k++) {
                inputs.H_s.xz(i, j, k) =
                        PSTD.ca(n, k) * inputs.H_s.xz(i, j, k) +
                        PSTD.cb(n, k) * eh_vec[n][k][0] / ((double) PSTD.N_hz);
              }
            }

          // PSTD, H_s.xz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else
         // PseudoSpectral)

        if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, H_s.xy
#pragma omp for
          // H_s.xy updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.xy(i, j, k) =
                          inputs.D.a.y[array_ind] * inputs.H_s.xy(i, j, k) +
                          inputs.D.b.y[array_ind] *
                                  (inputs.E_s.zy(i, j, k) +
                                   inputs.E_s.zx(i, j, k) -
                                   inputs.E_s.zy(i, j + 1, k) -
                                   inputs.E_s.zx(i, j + 1, k));
                else
                  inputs.H_s.xy(i, j, k) =
                          inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.xy(i, j, k) +
                          inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.zy(i, j, k) +
                                   inputs.E_s.zx(i, j, k) -
                                   inputs.E_s.zy(i, j + 1, k) -
                                   inputs.E_s.zx(i, j + 1, k));
              }
          // FDTD, H_s.xy
        } else {
#pragma omp for
          // H_s.xy updates
          for (k = 0; k < K_tot; k++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (j = 0; j < J_tot; j++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca(n, j) = inputs.D.a.y[array_ind];
                  PSTD.cb(n, j) = inputs.D.b.y[array_ind];
                } else {
                  PSTD.ca(n, j) =
                          inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1];
                  PSTD.cb(n, j) =
                          inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1];
                  //		H_s.xy(i,j,k) =
                  // Dmaterial.Da.y[materials[k][j][i]-1]*H_s.xy(i,j,k)+Dmaterial.Db.y[materials[k][j][i]-1]*(E_s.zy(i,j,k)
                  //+ E_s.zx(i,j,k) - E_s.zy[k][j+1][i] - E_s.zx[k][j+1][i]);
                }

                eh_vec[n][j][0] =
                        inputs.E_s.zy(i, j, k) + inputs.E_s.zx(i, j, k);
                eh_vec[n][j][1] = 0.;
              }
              j = J_tot;
              eh_vec[n][j][0] = inputs.E_s.zy(i, j, k) + inputs.E_s.zx(i, j, k);
              eh_vec[n][j][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hy, PSTD.N_hy,
                               inputs.H_s.xy.plan_f[n],
                               inputs.H_s.xy.plan_b[n]);

              for (j = 0; j < J_tot; j++) {
                inputs.H_s.xy(i, j, k) =
                        PSTD.ca(n, j) * inputs.H_s.xy(i, j, k) -
                        PSTD.cb(n, j) * eh_vec[n][j][0] / ((double) PSTD.N_hy);
              }
            }
          // PSTD, H_s.xy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else
         // PseudoSpectral)

        if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, H_s.yx
#pragma omp for
          // H_s.yx updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound_plus_1; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.yx(i, j, k) =
                          inputs.D.a.x[array_ind] * inputs.H_s.yx(i, j, k) +
                          inputs.D.b.x[array_ind] *
                                  (inputs.E_s.zx(i + 1, j, k) +
                                   inputs.E_s.zy(i + 1, j, k) -
                                   inputs.E_s.zx(i, j, k) -
                                   inputs.E_s.zy(i, j, k));
                else {
                  inputs.H_s.yx(i, j, k) =
                          inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.yx(i, j, k) +
                          inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.zx(i + 1, j, k) +
                                   inputs.E_s.zy(i + 1, j, k) -
                                   inputs.E_s.zx(i, j, k) -
                                   inputs.E_s.zy(i, j, k));
                }
              }
          // FDTD, H_s.yx
        } else {
#pragma omp for
          // H_s.yx updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound_plus_1; j++) {
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca(n, i) = inputs.D.a.x[array_ind];
                  PSTD.cb(n, i) = inputs.D.b.x[array_ind];
                } else {
                  PSTD.ca(n, i) =
                          inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1];
                  PSTD.cb(n, i) =
                          inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1];
                  //	H_s.yx(i, j, k) =
                  // Dmaterial.Da.x[materials[k][j][i]-1]*H_s.yx(i, j,
                  // k)+Dmaterial.Db.x[materials[k][j][i]-1]*(E_s.zx[k][j][i+1]
                  //+ E_s.zy[k][j][i+1] - E_s.zx(i,j,k) - E_s.zy(i,j,k));
                }

                eh_vec[n][i][0] =
                        inputs.E_s.zx(i, j, k) + inputs.E_s.zy(i, j, k);
                eh_vec[n][i][1] = 0.;
              }
              i = I_tot;
              eh_vec[n][i][0] = inputs.E_s.zx(i, j, k) + inputs.E_s.zy(i, j, k);
              eh_vec[n][i][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hx, PSTD.N_hx,
                               inputs.H_s.yx.plan_f[n],
                               inputs.H_s.yx.plan_b[n]);

              for (i = 0; i < I_tot; i++) {
                inputs.H_s.yx(i, j, k) =
                        PSTD.ca(n, i) * inputs.H_s.yx(i, j, k) +
                        PSTD.cb(n, i) * eh_vec[n][i][0] / ((double) PSTD.N_hx);
              }
            }
          // PSTD, H_s.yx
        }

        if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, H_s.yz
#pragma omp for
          // H_s.yz updates
          for (k = 0; k < K_tot; k++) {
            for (j = 0; j < loop_variables.J_loop_upper_bound_plus_1; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.materials[k][j][i]) {
                  inputs.H_s.yz(i, j, k) =
                          inputs.D.a.z[k_loc] * inputs.H_s.yz(i, j, k) +
                          inputs.D.b.z[k_loc] * (inputs.E_s.xy(i, j, k) +
                                                 inputs.E_s.xz(i, j, k) -
                                                 inputs.E_s.xy(i, j, k + 1) -
                                                 inputs.E_s.xz(i, j, k + 1));
                } else {
                  inputs.H_s.yz(i, j, k) =
                          inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.yz(i, j, k) +
                          inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.xy(i, j, k) +
                                   inputs.E_s.xz(i, j, k) -
                                   inputs.E_s.xy(i, j, k + 1) -
                                   inputs.E_s.xz(i, j, k + 1));
                }
              }
          }
          // FDTD, H_s.yz
        } else {
          // #pragma omp for
          // H_s.yz updates
          for (j = 0; j < loop_variables.J_loop_upper_bound_plus_1; j++)
#pragma omp for
            for (i = 0; i < I_tot; i++) {
              for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca(n, k) = inputs.D.a.z[k_loc];
                  PSTD.cb(n, k) = inputs.D.b.z[k_loc];
                } else {
                  PSTD.ca(n, k) =
                          inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1];
                  PSTD.cb(n, k) =
                          inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1];
                  // H_s.yz(i,j,k) =
                  // Dmaterial.Da.z[materials[k][j][i]-1]*H_s.yz(i,j,k)+Dmaterial.Db.z[materials[k][j][i]-1]*(E_s.xy(i,j,k)
                  // + E_s.xz(i, j, k) - E_s.xy[k+1][j][i] - E_s.xz[k+1][j][i]);
                }

                eh_vec[n][k][0] =
                        inputs.E_s.xy(i, j, k) + inputs.E_s.xz(i, j, k);
                eh_vec[n][k][1] = 0.;
              }
              k = K_tot;
              eh_vec[n][k][0] = inputs.E_s.xy(i, j, k) + inputs.E_s.xz(i, j, k);
              eh_vec[n][k][1] = 0.;
              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hz, PSTD.N_hz,
                               inputs.H_s.yz.plan_f[n],
                               inputs.H_s.yz.plan_b[n]);

              for (k = 0; k < K_tot; k++) {
                inputs.H_s.yz(i, j, k) =
                        PSTD.ca(n, k) * inputs.H_s.yz(i, j, k) -
                        PSTD.cb(n, k) * eh_vec[n][k][0] / ((double) PSTD.N_hz);
              }
            }
          // PSTD, H_s.yz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else
         // PseudoSpectral)
      }  //(params.dimension==THREE || params.dimension==TE)
      else {

#pragma omp for
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++)
              if (!inputs.materials[k][j][i]) inputs.H_s.xz(i, j, k) = 0.;
              else
                inputs.H_s.xz(i, j, k) = 0.;

#pragma omp for
        // H_s.xy update
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              k_loc = k;
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl &&
                    k < (inputs.params.pml.Dzl +
                         loop_variables.n_non_pml_cells_in_K)) {
                  if ((k - inputs.structure[i][1]) <
                              (loop_variables.n_non_pml_cells_in_K +
                               inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.n_non_pml_cells_in_K +
                            inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl +
                            loop_variables.n_non_pml_cells_in_K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = j;
              else
                array_ind = (J_tot + 1) * k_loc + j;
              if (!inputs.materials[k][j][i])
                inputs.H_s.xy(i, j, k) =
                        inputs.D.a.y[array_ind] * inputs.H_s.xy(i, j, k) +
                        inputs.D.b.y[array_ind] * (inputs.E_s.zy(i, j, k) +
                                                   inputs.E_s.zx(i, j, k) -
                                                   inputs.E_s.zy(i, j + 1, k) -
                                                   inputs.E_s.zx(i, j + 1, k));
              else
                inputs.H_s.xy(i, j, k) =
                        inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] *
                                inputs.H_s.xy(i, j, k) +
                        inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                (inputs.E_s.zy(i, j, k) +
                                 inputs.E_s.zx(i, j, k) -
                                 inputs.E_s.zy(i, j + 1, k) -
                                 inputs.E_s.zx(i, j + 1, k));
            }

#pragma omp for
        // H_s.yx update
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++) {
              k_loc = k;
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl &&
                    k < (inputs.params.pml.Dzl +
                         loop_variables.n_non_pml_cells_in_K)) {
                  if ((k - inputs.structure[i][1]) <
                              (loop_variables.n_non_pml_cells_in_K +
                               inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.n_non_pml_cells_in_K +
                            inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl +
                            loop_variables.n_non_pml_cells_in_K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = i;
              else
                array_ind = (I_tot + 1) * k_loc + i;
              if (!inputs.materials[k][j][i])
                inputs.H_s.yx(i, j, k) =
                        inputs.D.a.x[array_ind] * inputs.H_s.yx(i, j, k) +
                        inputs.D.b.x[array_ind] * (inputs.E_s.zx(i + 1, j, k) +
                                                   inputs.E_s.zy(i + 1, j, k) -
                                                   inputs.E_s.zx(i, j, k) -
                                                   inputs.E_s.zy(i, j, k));
              else
                inputs.H_s.yx(i, j, k) =
                        inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] *
                                inputs.H_s.yx(i, j, k) +
                        inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                (inputs.E_s.zx(i + 1, j, k) +
                                 inputs.E_s.zy(i + 1, j, k) -
                                 inputs.E_s.zx(i, j, k) -
                                 inputs.E_s.zy(i, j, k));
            }

#pragma omp for
        for (k = 0; k <= K_tot; k++) {
          for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++)
              if (!inputs.materials[k][j][i]) inputs.H_s.yz(i, j, k) = 0.;
              else
                inputs.H_s.yz(i, j, k) = 0.;
        }
      }

      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, H_s.zy
#pragma omp for
          // H_s.zy update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.zy(i, j, k) =
                          inputs.D.a.y[array_ind] * inputs.H_s.zy(i, j, k) +
                          inputs.D.b.y[array_ind] *
                                  (inputs.E_s.xy(i, j + 1, k) +
                                   inputs.E_s.xz(i, j + 1, k) -
                                   inputs.E_s.xy(i, j, k) -
                                   inputs.E_s.xz(i, j, k));
                else
                  inputs.H_s.zy(i, j, k) =
                          inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.zy(i, j, k) +
                          inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.xy(i, j + 1, k) +
                                   inputs.E_s.xz(i, j + 1, k) -
                                   inputs.E_s.xy(i, j, k) -
                                   inputs.E_s.xz(i, j, k));
              }
          // FDTD, H_s.zy
        } else {
#pragma omp for
          // H_s.zy update
          for (k = 0; k < (K_tot + 1); k++)
            for (i = 0; i < I_tot; i++) {
              for (j = 0; j < J_tot; j++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca(n, j) = inputs.D.a.y[array_ind];
                  PSTD.cb(n, j) = inputs.D.b.y[array_ind];
                } else {
                  PSTD.ca(n, j) =
                          inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1];
                  PSTD.cb(n, j) =
                          inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1];
                  //	      H_s.zy(i,j,k) =
                  // Dmaterial.Da.y[materials[k][j][i]-1]*H_s.zy(i,j,k)+Dmaterial.Db.y[materials[k][j][i]-1]*(E_s.xy[k][j+1][i]
                  //+ E_s.xz[k][j+1][i] - E_s.xy(i,j,k) - E_s.xz(i, j, k));
                }

                eh_vec[n][j][0] =
                        inputs.E_s.xy(i, j, k) + inputs.E_s.xz(i, j, k);
                eh_vec[n][j][1] = 0.;
              }
              j = J_tot;
              eh_vec[n][j][0] = inputs.E_s.xy(i, j, k) + inputs.E_s.xz(i, j, k);
              eh_vec[n][j][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hy, PSTD.N_hy,
                               inputs.H_s.zy.plan_f[n],
                               inputs.H_s.zy.plan_b[n]);

              for (j = 0; j < J_tot; j++) {
                inputs.H_s.zy(i, j, k) =
                        PSTD.ca(n, j) * inputs.H_s.zy(i, j, k) +
                        PSTD.cb(n, j) * eh_vec[n][j][0] / ((double) PSTD.N_hy);
              }
            }
          // PSTD, H_s.zy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else
         // PseudoSpectral)


        if (solver_method == SolverMethod::FiniteDifference) {
// FDTD, H_s.zx
#pragma omp for
          // H_s.zx update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.zx(i, j, k) =
                          inputs.D.a.x[array_ind] * inputs.H_s.zx(i, j, k) +
                          inputs.D.b.x[array_ind] *
                                  (inputs.E_s.yx(i, j, k) +
                                   inputs.E_s.yz(i, j, k) -
                                   inputs.E_s.yx(i + 1, j, k) -
                                   inputs.E_s.yz(i + 1, j, k));
                else
                  inputs.H_s.zx(i, j, k) =
                          inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.zx(i, j, k) +
                          inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.yx(i, j, k) +
                                   inputs.E_s.yz(i, j, k) -
                                   inputs.E_s.yx(i + 1, j, k) -
                                   inputs.E_s.yz(i + 1, j, k));
              }
          // FDTD, H_s.zx
        } else {
#pragma omp for
          // H_s.zx update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_loop_upper_bound; j++) {
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl &&
                      k < (inputs.params.pml.Dzl +
                           loop_variables.n_non_pml_cells_in_K)) {
                    if ((k - inputs.structure[i][1]) <
                                (loop_variables.n_non_pml_cells_in_K +
                                 inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.n_non_pml_cells_in_K +
                              inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl +
                              loop_variables.n_non_pml_cells_in_K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca(n, i) = inputs.D.a.x[array_ind];
                  PSTD.cb(n, i) = inputs.D.b.x[array_ind];
                } else {
                  PSTD.ca(n, i) =
                          inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1];
                  PSTD.cb(n, i) =
                          inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1];
                }

                eh_vec[n][i][0] =
                        inputs.E_s.yx(i, j, k) + inputs.E_s.yz(i, j, k);
                eh_vec[n][i][1] = 0.;
              }
              i = I_tot;
              eh_vec[n][i][0] = inputs.E_s.yx(i, j, k) + inputs.E_s.yz(i, j, k);
              eh_vec[n][i][1] = 0.;


              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hx, PSTD.N_hx,
                               inputs.H_s.zx.plan_f[n],
                               inputs.H_s.zx.plan_b[n]);

              for (i = 0; i < I_tot; i++) {
                inputs.H_s.zx(i, j, k) =
                        PSTD.ca(n, i) * inputs.H_s.zx(i, j, k) -
                        PSTD.cb(n, i) * eh_vec[n][i][0] / ((double) PSTD.N_hx);
              }
            }
          // PSTD, H_s.zx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else
         // PseudoSpectral)
      }  //(params.dimension==THREE || params.dimension==TE)
    }    // end parallel
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    /* Update source terms for self consistency across scattered/total interface
     * - H updates (use time_E) */

    if (inputs.params.source_mode == SourceMode::steadystate) {
      // Steady-state source term updates
      H_source_update_all_steadystate(time_E);

      //! Common phase term in update equations
      complex<double> common_phase =
              exp(-IMAGINARY_UNIT *
                  fmod(inputs.params.omega_an * time_E, 2. * DCPI));
      //! Common amplitude term in update equations
      double common_amplitude = linear_ramp(time_E);
      // Update output field
      outputs.E.ft = real(common_amplitude * common_phase);
    } else if (inputs.params.source_mode == SourceMode::pulsed) {
      // Pulsed source term updates
      update_source_terms_pulsed(time_E, tind);

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
      // Update output field
      outputs.E.ft = common_amplitude * real(common_phase);
    }
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    // If it is time for a new acquisition period, update the normalisation
    // factors and extract phasors
    new_acquisition_period(tind);
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    // Perform setup for next iteration
    end_of_iteration_steps(time_of_last_log_write, tind,
                           loop_variables.E_at_previous_iteration);
  }
  // end of main iteration loop

  if (TIME_MAIN_LOOP) {
    timers.end_timer(TimersTrackingLoop::MAIN);
    spdlog::info("Time elapsed in main loop (s): {0:e}",
                 timers.time_ellapsed_by(TimersTrackingLoop::MAIN));
  }
}
