#include "iterator_executor.h"

#include <complex>
#include <omp.h>
#include <spdlog/spdlog.h>

#include "globals.h"
#include "numerical_derivative.h"

using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;
using namespace tdms_phys_constants;
using namespace std;

void Iterator_Executor::run_main_loop() {
    // declare variables that we need whilst the main loop runs

    int n;      //< omp thread number (local to thread)
    int i, j, k;//< loop indices (local to omp threads when threaded)
    int k_loc;  //< indexing variable (local to thread)

    double rho;                     //< material constants (local to thread)
    double alpha_l, beta_l, gamma_l;//< material constants (local to thread)
    double kappa_l, sigma_l;        //< material constants (local to thread)

    // C and D vars for free-space and the perfectly matched layer (local to thread)
    double Enp1, Jnp1;
    double Ca, Cb, Cc;//< used by the interpolation scheme (local to thread)

    complex<double> Idxt, Idyt, kprop;

    double phaseTermE;               //< phase term (local to thread)
    std::complex<double> cphaseTermE;//< phase term (local to thread)
    double lambda_an_t;              //< wavelength in air (local to thread)

    if (TIME_MAIN_LOOP) { start_timer(IterationTimers::MAIN); }

    for (tind = params.start_tind; tind < params.Nt; tind++) {
        // Update the "time" the fields are currently at
        update_field_times_to_current_iteration(params.dt);

        //Extract phasors
        start_timer(IterationTimers::INTERNAL);
        if ((dft_counter == Nsteps) && (params.run_mode == RunMode::complete) &&
            (params.source_mode == SourceMode::steadystate) && params.exphasorsvolume) {

        dft_counter = 0;

        double tol = E.normalised_difference(E_copy);
        if (tol < TOL) break;//required accuracy obtained

        spdlog::debug("Phasor convergence: {} (actual) > {} (required)", tol, TOL);
        E_copy.set_values_from(E);

        E.zero();
        H.zero();
        spdlog::debug("Zeroed the phasors");

        if (params.exphasorssurface) {
            surface_phasors.zero_surface_EH();
            spdlog::debug("Zeroed the surface components");
        }
        }

        if ((params.source_mode == SourceMode::steadystate) &&
            (params.run_mode == RunMode::complete) && params.exphasorsvolume) {

        E.set_phasors(E_s, dft_counter - 1, params.omega_an, params.dt, Nsteps);
        H.set_phasors(H_s, dft_counter, params.omega_an, params.dt, Nsteps);

        if (params.exphasorssurface) {
            if (params.intphasorssurface) {
            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
                surface_phasors.extractPhasorsSurface(ifx, E_s, H_s, dft_counter,
                                                    f_ex_vec[ifx] * 2 * DCPI, Nsteps, params);
            }
            dft_counter++;
            } else {
            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
                // do not interpolate when extracting
                surface_phasors.extractPhasorsSurface(
                        ifx, E_s, H_s, dft_counter, f_ex_vec[ifx] * 2 * DCPI, Nsteps, params, false);
            }
            dft_counter++;
            }
        }

        } else if ((params.source_mode == SourceMode::pulsed) &&
                    (params.run_mode == RunMode::complete) && params.exphasorsvolume) {
        if (TIME_EXEC) {
            click_timer(IterationTimers::INTERNAL);
            ;
        }

        if ((tind - params.start_tind) % params.Np == 0) {
            E.set_phasors(E_s, tind - 1, params.omega_an, params.dt, params.Npe);
            H.set_phasors(H_s, tind, params.omega_an, params.dt, params.Npe);
        }
        if (TIME_EXEC) {
            click_timer(IterationTimers::INTERNAL);
            ;
        }
        //fprintf(stderr,"Pos 01b:\n");
        }

        /*extract fieldsample*/
        if (fieldsample.all_vectors_are_non_empty()) {
        fieldsample.extract(E_s, params.pml, params.Nt);
        }
        /*end extract fieldsample*/

        //fprintf(stderr,"Pos 02:\n");
        if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete &&
            params.exphasorssurface) {
        if ((tind - params.start_tind) % params.Np == 0) {
            if (params.intphasorssurface)
            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
                surface_phasors.extractPhasorsSurface(ifx, E_s, H_s, tind, f_ex_vec[ifx] * 2 * DCPI,
                                                    params.Npe, params);
            }
            else
            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
                // do not interpolate when extracting
                surface_phasors.extractPhasorsSurface(ifx, E_s, H_s, tind, f_ex_vec[ifx] * 2 * DCPI,
                                                    params.Npe, params, false);
            }
        }
        }

        if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete &&
            (vertex_phasors.there_are_vertices_to_extract_at()) &&
            ((tind - params.start_tind) % params.Np == 0)) {
        //     fprintf(stderr,"loc 01 (%d,%d,%d)\n",tind,params.start_tind,params.Np);
        //fprintf(stderr,"loc 03\n");
        //	  fprintf(stderr,"EPV 01\n");
        for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            vertex_phasors.extractPhasorsVertices(ifx, E_s, H_s, tind, f_ex_vec[ifx] * 2 * DCPI,
                                                params);
        }
        }

        //fprintf(stderr,"Pos 02a:\n");
        if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete &&
            params.exdetintegral) {
        if ((tind - params.start_tind) % params.Np == 0) {
            spdlog::debug("Setting Ex_t, Ey_t");

            //First need to sum up the Ex and Ey values on a plane ready for FFT, remember that Ex_t and
            // Ey_t are in row-major format whilst Exy etc. are in column major format
            for (j = params.pml.Dyl; j < (J_tot - params.pml.Dyu); j++)
            for (i = params.pml.Dxl; i < (I_tot - params.pml.Dxu); i++) {
                int m = j - params.pml.Dyl +
                        (i - params.pml.Dxl) * (J_tot - params.pml.Dyu - params.pml.Dyl);
                Ex_t.v[m][0] = E_s.xy[params.k_det_obs][j][i] + E_s.xz[params.k_det_obs][j][i];
                Ex_t.v[m][1] = 0.;
                Ey_t.v[m][0] = E_s.yx[params.k_det_obs][j][i] + E_s.yz[params.k_det_obs][j][i];
                Ey_t.v[m][1] = 0.;
            }

            fftw_execute(Ex_t.plan);
            fftw_execute(Ey_t.plan);

            //Iterate over each mode
            for (int im = 0; im < D_tilde.num_det_modes(); im++) {

            //Now go back to column-major
            for (j = 0; j < (J_tot - params.pml.Dyu - params.pml.Dyl); j++)
                for (i = 0; i < (I_tot - params.pml.Dxu - params.pml.Dxl); i++) {
                int m = j + i * (J_tot - params.pml.Dyu - params.pml.Dyl);
                Ex_t.cm[j][i] = Ex_t.v[m][0] + IMAGINARY_UNIT * Ex_t.v[m][1];
                Ey_t.cm[j][i] = Ey_t.v[m][0] + IMAGINARY_UNIT * Ey_t.v[m][1];
                }

            //fprintf(stderr,"Pos 02a [3]:\n");
            //Now multiply the pupil, mostly the pupil is non-zero in only a elements
            for (j = 0; j < (J_tot - params.pml.Dyu - params.pml.Dyl); j++)
                for (i = 0; i < (I_tot - params.pml.Dxu - params.pml.Dxl); i++) {
                Ex_t.cm[j][i] *= pupil[j][i] * D_tilde.x[j][i][im];
                Ey_t.cm[j][i] *= pupil[j][i] * D_tilde.y[j][i][im];
                }

                //now iterate over each frequency to extract phasors at
    #pragma omp parallel default(shared) private(lambda_an_t, Idxt, Idyt, i, j, kprop, phaseTermE,     \
                                                cphaseTermE)
            {
    #pragma omp for
            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
                //wavelength in air
                lambda_an_t = LIGHT_V / f_ex_vec[ifx];
                //fprintf(stdout,"lambda_an_t = %e, LIGHT_V = %e, z_obs = %e\n",lambda_an_t,LIGHT_V,z_obs);
                Idxt = 0.;
                Idyt = 0.;

                //now loop over all angular frequencies
                for (j = 0; j < (J_tot - params.pml.Dyu - params.pml.Dyl); j++)
                for (i = 0; i < (I_tot - params.pml.Dxu - params.pml.Dxl); i++) {
                    if ((lambda_an_t * f_vec.x[i] * lambda_an_t * f_vec.x[i] +
                        lambda_an_t * f_vec.y[j] * lambda_an_t * f_vec.y[j]) < 1) {

                    if (!params.air_interface_present) {
                        //This had to be fixed since we must take into account the refractive index of the medium.
                        kprop = exp(IMAGINARY_UNIT * params.z_obs * 2. * DCPI / lambda_an_t * refind *
                                    sqrt(1. - pow(lambda_an_t * f_vec.x[i] / refind, 2.) -
                                        pow(lambda_an_t * f_vec.y[j] / refind, 2.)));
                        //fprintf(stdout,"%d %d %e %e %e %e %e %e %e\n",i,j,f_vec.x[i],f_vec.y[j],real(kprop),imag(kprop),z_obs,DCPI,lambda_an_t);
                    } else {
                        kprop = exp(IMAGINARY_UNIT * (-params.air_interface + params.z_obs) * 2. *
                                    DCPI / lambda_an_t * refind *
                                    sqrt(1. - pow(lambda_an_t * f_vec.x[i] / refind, 2.) -
                                        pow(lambda_an_t * f_vec.y[j] / refind, 2.))) *
                                exp(IMAGINARY_UNIT * params.air_interface * 2. * DCPI / lambda_an_t *
                                    sqrt(1. - pow(lambda_an_t * f_vec.x[i], 2.) -
                                        pow(lambda_an_t * f_vec.y[j], 2.)));
                    }
                    } else
                    kprop = 0.;

                    Idxt += Ex_t.cm[j][i] * kprop;
                    Idyt += Ey_t.cm[j][i] * kprop;
                }
                phaseTermE = fmod(f_ex_vec[ifx] * 2. * DCPI * ((double) tind) * params.dt, 2 * DCPI);
                cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. / ((double) params.Npe);

                Idx[ifx][im] += Idxt * cphaseTermE;
                Idy[ifx][im] += Idyt * cphaseTermE;

            }//end of loop on frequencies
            }  //end of pragma omp parallel
        }    //end of loop over each mode
        }
    }//end of section for calculating detector function

    //fprintf(stderr,"Pos 02b:\n");
    if (params.run_mode == RunMode::complete)
        if (params.dimension == THREE) {
        // extract the phasors just above the line
        extract_phasors_in_plane();
        }
    //fprintf(stderr,"Pos 02c:\n");

    //Update equations for the E field

    /*There are two options for determining the update coefficients for the FDTD cell:

        1) If cell (i,j,k) is either free space or PML:

        materials[k][j][i] will be set to 0. In this case the update parameter used will
        be given by C.a.y[j], C.b.y[j] etc depending on which update equation is being implemented.

        2) if cell (i,j,k) is composed of a scattering type material then materials[k][j][i] will be
        non-zero and will be an index into Cmaterial.a.y and Cmaterial.b.y etc depending on which
        update equation is being implemented.

    */

    int array_ind = 0;
    //fprintf(stderr,"I_tot=%d, J_tot=%d, K_tot=%d\n",I_tot,J_tot,K_tot);
    if (TIME_EXEC) {
        click_timer(IterationTimers::INTERNAL);
        ;
    }
    //fprintf(stderr,"Dimension = %d\n",params.dimension);
    /*
        for(k=0;k<(K_tot+1);k++)
        fprintf(stdout,"%e ",Exy[k][13][13]+Exz[k][13][13]);
        fprintf(stdout,"\n");
    */
    (void) n;// n is unused in FD derivatives – this silences the compiler warning

    #pragma omp parallel default(shared) private(i, j, k, n, rho, k_loc, array_ind, Ca, Cb, Cc,        \
                                                alpha_l, beta_l, gamma_l, kappa_l, sigma_l, Enp1,     \
                                                Jnp1)//,ca_vec,cb_vec,cc_vec,eh_vec)
    {
        n = omp_get_thread_num();
        Enp1 = 0.0;
        array_ind = 0;

        if (params.dimension == THREE || params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
            //FDTD, E_s.xy
    #pragma omp for
            for (k = 0; k < (K_tot + 1); k++)
            for (j = 1; j < J_tot; j++)
                for (i = 0; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                    //fprintf(stdout,"(%d,%d,%d,%d)\n",i,j,k,tind);
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.y[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.y[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.y[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k][j][i + 1]) {
                        Ca = Ca + C.a.y[array_ind];
                        Cb = Cb + C.b.y[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.y[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.y[materials[k][j][i + 1] - 1];
                        Cb = Cb + Cmaterial.b.y[materials[k][j][i + 1] - 1];
                        Cc = Cc + Cmaterial.c.y[materials[k][j][i + 1] - 1];
                    }
                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.y[array_ind];
                    kappa_l = matched_layer.kappa.y[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][j][i + 1]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][j][i + 1]) {
                        alpha_l += alpha[materials[k][j][i + 1] - 1];
                        beta_l += beta[materials[k][j][i + 1] - 1];
                        gamma_l += gamma[materials[k][j][i + 1] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }


                Enp1 = Ca * E_s.xy[k][j][i] + Cb * (H_s.zy[k][j][i] + H_s.zx[k][j][i] -
                                                    H_s.zy[k][j - 1][i] - H_s.zx[k][j - 1][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.xy[k][j][i] -
                            1. / 2. * Cb * params.delta.dy *
                                    ((1 + alpha_l) * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dy * J_c.xy[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xy[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.xy[k][j][i];

                    E_nm1.xy[k][j][i] = E_s.xy[k][j][i];
                    J_nm1.xy[k][j][i] = J_s.xy[k][j][i];
                    J_s.xy[k][j][i] = Jnp1;

                    //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,J_s.xy[k][j][i]);
                }

                if (is_conductive && rho) { J_c.xy[k][j][i] -= rho * (Enp1 + E_s.xy[k][j][i]); }

                E_s.xy[k][j][i] = Enp1;
                }
            //FDTD, E_s.xy
        } else {
    //fprintf(stderr,"Pos 02d:\n");
    #pragma omp for
            for (k = 0; k < (K_tot + 1); k++)
            for (i = 0; i < I_tot; i++) {
                for (j = 1; j < J_tot; j++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                    //fprintf(stdout,"(%d,%d,%d,%d)\n",i,j,k,tind);
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.y[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.y[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.y[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k][j][i + 1]) {
                        Ca = Ca + C.a.y[array_ind];
                        Cb = Cb + C.b.y[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.y[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.y[materials[k][j][i + 1] - 1];
                        Cb = Cb + Cmaterial.b.y[materials[k][j][i + 1] - 1];
                        Cc = Cc + Cmaterial.c.y[materials[k][j][i + 1] - 1];
                    }
                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.y[array_ind];
                    kappa_l = matched_layer.kappa.y[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][j][i + 1]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][j][i + 1]) {
                        alpha_l += alpha[materials[k][j][i + 1] - 1];
                        beta_l += beta[materials[k][j][i + 1] - 1];
                        gamma_l += gamma[materials[k][j][i + 1] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }


                Enp1 = 0.0;
                //Enp1 = Ca*E_s.xy[k][j][i]+Cb*(H_s.zy[k][j][i] + H_s.zx[k][j][i] - H_s.zy[k][j-1][i] - H_s.zx[k][j-1][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.xy[k][j][i] -
                            1. / 2. * Cb * params.delta.dy *
                                    ((1 + alpha_l) * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dy * J_c.xy[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xy[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.xy[k][j][i];

                    E_nm1.xy[k][j][i] = E_s.xy[k][j][i];
                    J_nm1.xy[k][j][i] = J_s.xy[k][j][i];
                    J_s.xy[k][j][i] = Jnp1;

                    //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,J_s.xy[k][j][i]);
                }

                if (is_conductive && rho) { J_c.xy[k][j][i] -= rho * (Enp1 + E_s.xy[k][j][i]); }

                eh_vec[n][j][0] = H_s.zy[k][j][i] + H_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;
                ca_vec[n][j - 1] = Ca;
                cb_vec[n][j - 1] = Cb;
                }
                if (J_tot > 1) {
                j = 0;
                eh_vec[n][j][0] = H_s.zy[k][j][i] + H_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;
                first_derivative(eh_vec[n], eh_vec[n], dk_e_y, N_e_y, E_s.xy.plan_f[n],
                                    E_s.xy.plan_b[n]);


                //fprintf(stdout,"(%d,%d) %d (of %d)\n",i,k,n,omp_get_num_threads());

                for (j = 1; j < J_tot; j++) {
                    E_s.xy[k][j][i] = ca_vec[n][j - 1] * E_s.xy[k][j][i] +
                                    cb_vec[n][j - 1] * eh_vec[n][j][0] / ((double) N_e_y);
                }
                }
            }
            //PSTD, E_s.xy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
        /*
    if(is_disp){
    i=36;
    j=36;
    k=36;

    fprintf(stdout,"%e %e",J_s.xy[k][j][i],E_s.xy[k][j][i]);
    }
    */

        //fprintf(stderr,"Pos 04:\n");
        //E_s.xz updates
        if (solver_method == SolverMethod::FiniteDifference) {
    #pragma omp for
            for (k = 1; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++)
                for (i = 0; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.z[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.z[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.z[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k][j][i + 1]) {
                        Ca = Ca + C.a.z[k_loc];
                        Cb = Cb + C.b.z[k_loc];
                        if (params.is_disp_ml) Cc = Cc + C.c.z[k_loc];
                    } else {
                        Ca = Ca + Cmaterial.a.z[materials[k][j][i + 1] - 1];
                        Cb = Cb + Cmaterial.b.z[materials[k][j][i + 1] - 1];
                        Cc = Cc + Cmaterial.c.z[materials[k][j][i + 1] - 1];
                    }
                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.z[k_loc];
                    kappa_l = matched_layer.kappa.z[k_loc];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][j][i + 1]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][j][i + 1]) {
                        alpha_l += alpha[materials[k][j][i + 1] - 1];
                        beta_l += beta[materials[k][j][i + 1] - 1];
                        gamma_l += gamma[materials[k][j][i + 1] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }
                /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);
        if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Ca, Cb);*/
                Enp1 = Ca * E_s.xz[k][j][i] + Cb * (H_s.yx[k - 1][j][i] + H_s.yz[k - 1][j][i] -
                                                    H_s.yx[k][j][i] - H_s.yz[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.xz[k][j][i] -
                            1. / 2. * Cb * params.delta.dz *
                                    ((1 + alpha_l) * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dz * J_c.xz[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xz[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.xz[k][j][i];
                    E_nm1.xz[k][j][i] = E_s.xz[k][j][i];
                    J_nm1.xz[k][j][i] = J_s.xz[k][j][i];
                    J_s.xz[k][j][i] = Jnp1;
                }

                if (is_conductive && rho) { J_c.xz[k][j][i] -= rho * (Enp1 + E_s.xz[k][j][i]); }

                E_s.xz[k][j][i] = Enp1;
                }
            //FDTD, E_s.xz
        } else {
            //#pragma omp for
            for (j = 0; j < J_tot_p1_bound; j++)
    #pragma omp for
            for (i = 0; i < I_tot; i++) {
                for (k = 1; k < K_tot; k++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.z[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.z[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.z[materials[k][j][i] - 1];
                    }
                    if (params.interp_mat_props) {
                    if (!materials[k][j][i + 1]) {
                        Ca = Ca + C.a.z[k_loc];
                        Cb = Cb + C.b.z[k_loc];
                        if (params.is_disp_ml) Cc = Cc + C.c.z[k_loc];
                    } else {
                        Ca = Ca + Cmaterial.a.z[materials[k][j][i + 1] - 1];
                        Cb = Cb + Cmaterial.b.z[materials[k][j][i + 1] - 1];
                        Cc = Cc + Cmaterial.c.z[materials[k][j][i + 1] - 1];
                    }
                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.z[k_loc];
                    kappa_l = matched_layer.kappa.z[k_loc];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][j][i + 1]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][j][i + 1]) {
                        alpha_l += alpha[materials[k][j][i + 1] - 1];
                        beta_l += beta[materials[k][j][i + 1] - 1];
                        gamma_l += gamma[materials[k][j][i + 1] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }
                /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);
        if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Ca, Cb);*/
                //Enp1 = Ca*E_s.xz[k][j][i]+Cb*(H_s.yx[k-1][j][i] + H_s.yz[k-1][j][i] - H_s.yx[k][j][i] - H_s.yz[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.xz[k][j][i] -
                            1. / 2. * Cb * params.delta.dz *
                                    ((1 + alpha_l) * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dz * J_c.xz[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xz[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.xz[k][j][i];
                    E_nm1.xz[k][j][i] = E_s.xz[k][j][i];
                    J_nm1.xz[k][j][i] = J_s.xz[k][j][i];
                    J_s.xz[k][j][i] = Jnp1;
                }

                if (is_conductive && rho) { J_c.xz[k][j][i] -= rho * (Enp1 + E_s.xz[k][j][i]); }

                eh_vec[n][k][0] = H_s.yx[k][j][i] + H_s.yz[k][j][i];
                eh_vec[n][k][1] = 0.;
                ca_vec[n][k - 1] = Ca;
                cb_vec[n][k - 1] = Cb;
                }
                k = 0;
                eh_vec[n][k][0] = H_s.yx[k][j][i] + H_s.yz[k][j][i];
                eh_vec[n][k][1] = 0.;
                /*
    if (tind==1 & i==25 & j==25){
    for(k=0;k<N_e_z;k++)
    fprintf(stdout,"%e ",dk_e_z[k][0]);
    fprintf(stdout,"\n\n");
    for(k=0;k<N_e_z;k++)
    fprintf(stdout,"%e ",dk_e_z[k][1]);
    fprintf(stdout,"\n\n");
    for(k=0;k<N_e_z;k++)
    fprintf(stdout,"%e ",eh_vec[n][k][0]);
    fprintf(stdout,"\n\n");
    for(k=0;k<N_e_z;k++)
    fprintf(stdout,"%e ",eh_vec[n][k][1]);
    }
        */
                first_derivative(eh_vec[n], eh_vec[n], dk_e_z, N_e_z, E_s.xz.plan_f[n],
                                E_s.xz.plan_b[n]);
                /*
    if (tind==1 & i==25 & j==25){
    fprintf(stdout,"\n\n");
    for(k=0;k<N_e_z;k++)
    fprintf(stdout,"%e ",eh_vec[k][0]);
    fprintf(stdout,"\n\n");
    for(k=0;k<N_e_z;k++)
    fprintf(stdout,"%e ",eh_vec[k][1]);
    }
        */
                for (k = 1; k < K_tot; k++) {
                E_s.xz[k][j][i] = ca_vec[n][k - 1] * E_s.xz[k][j][i] -
                                    cb_vec[n][k - 1] * eh_vec[n][k][0] / ((double) N_e_z);
                }
            }
            //PSTD, E_s.xz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)

        //fprintf(stderr,"Pos 05:\n");
        //E_s.yx updates
        if (solver_method == SolverMethod::FiniteDifference) {
            //FDTD, E_s.yx
    #pragma omp for
            for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot_bound; j++)
                for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure) {
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                        Cc = 0;
                    } else {
                    Ca = Cmaterial.a.x[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.x[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.x[materials[k][j][i] - 1];
                    }
                    if (params.interp_mat_props) {
                    if (!materials[k][min(J_tot, j + 1)][i]) {
                        Ca = Ca + C.a.x[array_ind];
                        Cb = Cb + C.b.x[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.x[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.x[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cb = Cb + Cmaterial.b.x[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cc = Cc + Cmaterial.c.x[materials[k][min(J_tot, j + 1)][i] - 1];
                    }

                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.x[array_ind];
                    kappa_l = matched_layer.kappa.x[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][min(J_tot, j + 1)][i]) {
                        alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                        beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                        gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }


                Enp1 = Ca * E_s.yx[k][j][i] + Cb * (H_s.zx[k][j][i - 1] + H_s.zy[k][j][i - 1] -
                                                    H_s.zx[k][j][i] - H_s.zy[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.yx[k][j][i] -
                            1. / 2. * Cb * params.delta.dx *
                                    ((1 + alpha_l) * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dx * J_c.yx[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yx[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.yx[k][j][i];
                    E_nm1.yx[k][j][i] = E_s.yx[k][j][i];
                    J_nm1.yx[k][j][i] = J_s.yx[k][j][i];
                    J_s.yx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yx[k][j][i] -= rho * (Enp1 + E_s.yx[k][j][i]); }

                E_s.yx[k][j][i] = Enp1;
                }
            //FDTD, E_s.yx
        } else {
    #pragma omp for
            for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot_bound; j++) {
                for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure) {
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                        Cc = 0;
                    } else {
                    Ca = Cmaterial.a.x[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.x[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.x[materials[k][j][i] - 1];
                    }
                    if (params.interp_mat_props) {
                    if (!materials[k][min(J_tot, j + 1)][i]) {
                        Ca = Ca + C.a.x[array_ind];
                        Cb = Cb + C.b.x[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.x[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.x[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cb = Cb + Cmaterial.b.x[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cc = Cc + Cmaterial.c.x[materials[k][min(J_tot, j + 1)][i] - 1];
                    }

                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.x[array_ind];
                    kappa_l = matched_layer.kappa.x[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][min(J_tot, j + 1)][i]) {
                        alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                        beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                        gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }


                //Enp1 = Ca*E_s.yx[k][j][i]+Cb*(H_s.zx[k][j][i-1] + H_s.zy[k][j][i-1] - H_s.zx[k][j][i] - H_s.zy[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.yx[k][j][i] -
                            1. / 2. * Cb * params.delta.dx *
                                    ((1 + alpha_l) * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dx * J_c.yx[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yx[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.yx[k][j][i];
                    E_nm1.yx[k][j][i] = E_s.yx[k][j][i];
                    J_nm1.yx[k][j][i] = J_s.yx[k][j][i];
                    J_s.yx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yx[k][j][i] -= rho * (Enp1 + E_s.yx[k][j][i]); }

                eh_vec[n][i][0] = H_s.zx[k][j][i] + H_s.zy[k][j][i];
                eh_vec[n][i][1] = 0.;
                ca_vec[n][i - 1] = Ca;
                cb_vec[n][i - 1] = Cb;
                }
                i = 0;
                eh_vec[n][i][0] = H_s.zx[k][j][i] + H_s.zy[k][j][i];
                eh_vec[n][i][1] = 0.;

                first_derivative(eh_vec[n], eh_vec[n], dk_e_x, N_e_x, E_s.yx.plan_f[n],
                                E_s.yx.plan_b[n]);

                for (i = 1; i < I_tot; i++) {
                E_s.yx[k][j][i] = ca_vec[n][i - 1] * E_s.yx[k][j][i] -
                                    cb_vec[n][i - 1] * eh_vec[n][i][0] / ((double) N_e_x);
                //E_s.yx[k][j][i] = Enp1;
                }
            }
            //PSTD, E_s.yx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)

        //fprintf(stderr,"Pos 06:\n");
        //E_s.yz updates
        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, E_s.yz
    #pragma omp for
            for (k = 1; k < K_tot; k++)
            for (j = 0; j < J_tot_bound; j++)
                for (i = 0; i < (I_tot + 1); i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.z[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.z[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.z[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k][min(J_tot, j + 1)][i]) {
                        Ca = Ca + C.a.z[k_loc];
                        Cb = Cb + C.b.z[k_loc];
                        if (params.is_disp_ml) Cc = Cc + C.c.z[k_loc];
                    } else {
                        Ca = Ca + Cmaterial.a.z[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cb = Cb + Cmaterial.b.z[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cc = Cc + Cmaterial.c.z[materials[k][min(J_tot, j + 1)][i] - 1];
                    }

                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.z[k_loc];
                    kappa_l = matched_layer.kappa.z[k_loc];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][min(J_tot, j + 1)][i]) {
                        alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                        beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                        gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }

                //fprintf(stderr,"[%d %d %d]Ca: %e, Cb: %e, Cc: %e, alpha: %e, beta: %e, gamme: %e\n",i,j,k,Ca,Cb,Cc,alpha_l,beta_l,gamma_l);
                Enp1 = Ca * E_s.yz[k][j][i] + Cb * (H_s.xy[k][j][i] + H_s.xz[k][j][i] -
                                                    H_s.xy[k - 1][j][i] - H_s.xz[k - 1][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.yz[k][j][i] -
                            1. / 2. * Cb * params.delta.dz *
                                    ((1 + alpha_l) * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dz * J_c.yz[k][j][i];

                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yz[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.yz[k][j][i];
                    E_nm1.yz[k][j][i] = E_s.yz[k][j][i];
                    J_nm1.yz[k][j][i] = J_s.yz[k][j][i];
                    J_s.yz[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yz[k][j][i] -= rho * (Enp1 + E_s.yz[k][j][i]); }

                E_s.yz[k][j][i] = Enp1;
                }
            //FDTD, E_s.yz
        } else {
    #pragma omp for
            for (j = 0; j < J_tot_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
                for (k = 1; k < K_tot; k++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.z[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.z[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.z[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k][min(J_tot, j + 1)][i]) {
                        Ca = Ca + C.a.z[k_loc];
                        Cb = Cb + C.b.z[k_loc];
                        if (params.is_disp_ml) Cc = Cc + C.c.z[k_loc];
                    } else {
                        Ca = Ca + Cmaterial.a.z[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cb = Cb + Cmaterial.b.z[materials[k][min(J_tot, j + 1)][i] - 1];
                        Cc = Cc + Cmaterial.c.z[materials[k][min(J_tot, j + 1)][i] - 1];
                    }

                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.z[k_loc];
                    Cb = C.b.z[k_loc];
                    if (params.is_disp_ml) Cc = C.c.z[k_loc];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.z[k_loc];
                    kappa_l = matched_layer.kappa.z[k_loc];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k][min(J_tot, j + 1)][i]) {
                        alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                        beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                        gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }

                //fprintf(stderr,"[%d %d %d]Ca: %e, Cb: %e, Cc: %e, alpha: %e, beta: %e, gamme: %e\n",i,j,k,Ca,Cb,Cc,alpha_l,beta_l,gamma_l);
                //Enp1 = Ca*E_s.yz[k][j][i]+Cb*(H_s.xy[k][j][i] + H_s.xz[k][j][i] - H_s.xy[k-1][j][i] - H_s.xz[k-1][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.yz[k][j][i] -
                            1. / 2. * Cb * params.delta.dz *
                                    ((1 + alpha_l) * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dz * J_c.yz[k][j][i];

                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yz[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.yz[k][j][i];
                    E_nm1.yz[k][j][i] = E_s.yz[k][j][i];
                    J_nm1.yz[k][j][i] = J_s.yz[k][j][i];
                    J_s.yz[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yz[k][j][i] -= rho * (Enp1 + E_s.yz[k][j][i]); }

                eh_vec[n][k][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
                eh_vec[n][k][1] = 0.;
                ca_vec[n][k - 1] = Ca;
                cb_vec[n][k - 1] = Cb;
                }
                k = 0;
                eh_vec[n][k][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
                eh_vec[n][k][1] = 0.;
                first_derivative(eh_vec[n], eh_vec[n], dk_e_z, N_e_z, E_s.yz.plan_f[n],
                                E_s.yz.plan_b[n]);


                for (k = 1; k < K_tot; k++) {
                E_s.yz[k][j][i] = ca_vec[n][k - 1] * E_s.yz[k][j][i] +
                                    cb_vec[n][k - 1] * eh_vec[n][k][0] / ((double) N_e_z);
                //E_s.yz[k][j][i] = Enp1;
                }
            }
            //PSTD, E_s.yz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
        }  //if(params.dimension==THREE || params.dimension==TE)

        //fprintf(stderr,"Pos 07:\n");
        if (params.dimension == THREE || params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
    #pragma omp for
            //E_s.zx updates
            for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++)
                for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.x[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.x[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.x[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k + 1][j][i]) {
                        Ca = Ca + C.a.x[array_ind];
                        Cb = Cb + C.b.x[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.x[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.x[materials[k + 1][j][i] - 1];
                        Cb = Cb + Cmaterial.b.x[materials[k + 1][j][i] - 1];
                        Cc = Cc + Cmaterial.c.x[materials[k + 1][j][i] - 1];
                    }

                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.x[array_ind];
                    kappa_l = matched_layer.kappa.x[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k + 1][j][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k + 1][j][i]) {
                        alpha_l += alpha[materials[k + 1][j][i] - 1];
                        beta_l += beta[materials[k + 1][j][i] - 1];
                        gamma_l += gamma[materials[k + 1][j][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }

                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }

                /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);*/
                Enp1 = Ca * E_s.zx[k][j][i] + Cb * (H_s.yx[k][j][i] + H_s.yz[k][j][i] -
                                                    H_s.yx[k][j][i - 1] - H_s.yz[k][j][i - 1]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.zx[k][j][i] -
                            1. / 2. * Cb * params.delta.dx *
                                    ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dx * J_c.zx[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.zx[k][j][i];
                    E_nm1.zx[k][j][i] = E_s.zx[k][j][i];
                    J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                    J_s.zx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + E_s.zx[k][j][i]); }

                E_s.zx[k][j][i] = Enp1;
                }
            //FDTD, E_s.zx
        } else {
    #pragma omp for
            for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++) {
                for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.x[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.x[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.x[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k + 1][j][i]) {
                        Ca = Ca + C.a.x[array_ind];
                        Cb = Cb + C.b.x[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.x[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.x[materials[k + 1][j][i] - 1];
                        Cb = Cb + Cmaterial.b.x[materials[k + 1][j][i] - 1];
                        Cc = Cc + Cmaterial.c.x[materials[k + 1][j][i] - 1];
                    }

                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }
                } else {
                    Ca = C.a.x[array_ind];
                    Cb = C.b.x[array_ind];
                    if (params.is_disp_ml) Cc = C.c.x[array_ind];
                    else
                    Cc = 0.;
                    if (is_conductive) rho = rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.x[array_ind];
                    kappa_l = matched_layer.kappa.x[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k + 1][j][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k + 1][j][i]) {
                        alpha_l += alpha[materials[k + 1][j][i] - 1];
                        beta_l += beta[materials[k + 1][j][i] - 1];
                        gamma_l += gamma[materials[k + 1][j][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }

                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }

                /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);*/
                //Enp1 = Ca*E_s.zx[k][j][i]+Cb*(H_s.yx[k][j][i] + H_s.yz[k][j][i] - H_s.yx[k][j][i-1] - H_s.yz[k][j][i-1]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.zx[k][j][i] -
                            1. / 2. * Cb * params.delta.dx *
                                    ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dx * J_c.zx[k][j][i];
                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.zx[k][j][i];
                    E_nm1.zx[k][j][i] = E_s.zx[k][j][i];
                    J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                    J_s.zx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + E_s.zx[k][j][i]); }

                eh_vec[n][i][0] = H_s.yx[k][j][i] + H_s.yz[k][j][i];
                eh_vec[n][i][1] = 0.;
                ca_vec[n][i - 1] = Ca;
                cb_vec[n][i - 1] = Cb;
                }
                i = 0;
                eh_vec[n][i][0] = H_s.yx[k][j][i] + H_s.yz[k][j][i];
                eh_vec[n][i][1] = 0.;

                first_derivative(eh_vec[n], eh_vec[n], dk_e_x, N_e_x, E_s.zx.plan_f[n],
                                E_s.zx.plan_b[n]);

                for (i = 1; i < I_tot; i++) {
                E_s.zx[k][j][i] = ca_vec[n][i - 1] * E_s.zx[k][j][i] +
                                    cb_vec[n][i - 1] * eh_vec[n][i][0] / ((double) N_e_x);
                //E_s.zx[k][j][i] = Enp1;
                }
            }
            //PSTD, E_s.zx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
        }  //(params.dimension==THREE || params.dimension==TE)
        else {
    #pragma omp for
        //E_s.zx updates
        for (k = 0; k <= K_tot; k++)
            for (j = 0; j < (J_tot + 1); j++)
            for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                    k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                    k_loc = params.pml.Dzl + K - 1;
                    else
                    k_loc = params.pml.Dzl + 1;
                }
                if (!params.is_multilayer) array_ind = i;
                else
                array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (!materials[k][j][i]) {
                Ca = C.a.x[array_ind];
                Cb = C.b.x[array_ind];
                if (params.is_disp_ml) Cc = C.c.x[array_ind];
                else
                    Cc = 0.;
                if (is_conductive) rho = rho_cond.x[i];
                } else {
                rho = 0.;
                Ca = Cmaterial.a.x[materials[k][j][i] - 1];
                Cb = Cmaterial.b.x[materials[k][j][i] - 1];
                Cc = Cmaterial.c.x[materials[k][j][i] - 1];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;


                if (is_disp || params.is_disp_ml) {
                sigma_l = matched_layer.sigma.x[array_ind];
                kappa_l = matched_layer.kappa.x[array_ind];
                alpha_l = matched_layer.alpha[k_loc];
                beta_l = matched_layer.beta[k_loc];
                gamma_l = matched_layer.gamma[k_loc];

                if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];

                } else {
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                }
                }

                Enp1 = Ca * E_s.zx[k][j][i] + Cb * (H_s.yx[k][j][i] + H_s.yz[k][j][i] -
                                                    H_s.yx[k][j][i - 1] - H_s.yz[k][j][i - 1]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zx[k][j][i] -
                        1. / 2. * Cb * params.delta.dx *
                                ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dx * J_c.zx[k][j][i];

                if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                        kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.zx[k][j][i];
                E_nm1.zx[k][j][i] = E_s.zx[k][j][i];
                J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                J_s.zx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + E_s.zx[k][j][i]); }

                E_s.zx[k][j][i] = Enp1;
            }
        }
        //fprintf(stderr,"Pos 08:\n");
        if (params.dimension == THREE || params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
            //FDTD, E_s.zy
    #pragma omp for
            //E_s.zy updates
            for (k = 0; k < K_tot; k++)
            for (j = 1; j < J_tot; j++)
                for (i = 0; i < (I_tot + 1); i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.y[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.y[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.y[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k + 1][j][i]) {
                        Ca = Ca + C.a.y[array_ind];
                        Cb = Cb + C.b.y[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.y[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.y[materials[k + 1][j][i] - 1];
                        Cb = Cb + Cmaterial.b.y[materials[k + 1][j][i] - 1];
                        Cc = Cc + Cmaterial.c.y[materials[k + 1][j][i] - 1];
                    }
                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }

                } else {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                    Cc = 0;
                    if (is_conductive) rho = rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.y[array_ind];
                    kappa_l = matched_layer.kappa.y[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k + 1][j][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k + 1][j][i]) {
                        alpha_l += alpha[materials[k + 1][j][i] - 1];
                        beta_l += beta[materials[k + 1][j][i] - 1];
                        gamma_l += gamma[materials[k + 1][j][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }


                Enp1 = Ca * E_s.zy[k][j][i] + Cb * (H_s.xy[k][j - 1][i] + H_s.xz[k][j - 1][i] -
                                                    H_s.xy[k][j][i] - H_s.xz[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.zy[k][j][i] -
                            1. / 2. * Cb * params.delta.dy *
                                    ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dy * J_c.zy[k][j][i];

                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.zy[k][j][i];
                    E_nm1.zy[k][j][i] = E_s.zy[k][j][i];
                    J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                    J_s.zy[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + E_s.zy[k][j][i]); }
                E_s.zy[k][j][i] = Enp1;
                }
            //FDTD, E_s.zy
        } else {
    #pragma omp for
            //E_s.zy updates
            for (k = 0; k < K_tot; k++)
            for (i = 0; i < (I_tot + 1); i++) {
                for (j = 1; j < J_tot; j++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                    rho = 0.;
                    if (!materials[k][j][i]) {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                        Cc = 0.;
                    } else {
                    Ca = Cmaterial.a.y[materials[k][j][i] - 1];
                    Cb = Cmaterial.b.y[materials[k][j][i] - 1];
                    Cc = Cmaterial.c.y[materials[k][j][i] - 1];
                    }

                    if (params.interp_mat_props) {
                    if (!materials[k + 1][j][i]) {
                        Ca = Ca + C.a.y[array_ind];
                        Cb = Cb + C.b.y[array_ind];
                        if (params.is_disp_ml) Cc = Cc + C.c.y[array_ind];
                    } else {
                        Ca = Ca + Cmaterial.a.y[materials[k + 1][j][i] - 1];
                        Cb = Cb + Cmaterial.b.y[materials[k + 1][j][i] - 1];
                        Cc = Cc + Cmaterial.c.y[materials[k + 1][j][i] - 1];
                    }
                    Ca = Ca / 2.;
                    Cb = Cb / 2.;
                    Cc = Cc / 2.;
                    }

                } else {
                    Ca = C.a.y[array_ind];
                    Cb = C.b.y[array_ind];
                    if (params.is_disp_ml) Cc = C.c.y[array_ind];
                    else
                    Cc = 0;
                    if (is_conductive) rho = rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                    sigma_l = matched_layer.sigma.y[array_ind];
                    kappa_l = matched_layer.kappa.y[array_ind];
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                    if (materials[k][j][i] || materials[k + 1][j][i]) {
                    if (materials[k][j][i]) {
                        alpha_l = alpha[materials[k][j][i] - 1];
                        beta_l = beta[materials[k][j][i] - 1];
                        gamma_l = gamma[materials[k][j][i] - 1];
                    } else {
                        alpha_l = matched_layer.alpha[k_loc];
                        beta_l = matched_layer.beta[k_loc];
                        gamma_l = matched_layer.gamma[k_loc];
                    }

                    if (materials[k + 1][j][i]) {
                        alpha_l += alpha[materials[k + 1][j][i] - 1];
                        beta_l += beta[materials[k + 1][j][i] - 1];
                        gamma_l += gamma[materials[k + 1][j][i] - 1];
                    } else {
                        alpha_l += matched_layer.alpha[k_loc];
                        beta_l += matched_layer.beta[k_loc];
                        gamma_l += matched_layer.gamma[k_loc];
                    }
                    alpha_l = alpha_l / 2.;
                    beta_l = beta_l / 2.;
                    gamma_l = gamma_l / 2.;
                    }
                }


                //Enp1 = Ca*E_s.zy[k][j][i]+Cb*(H_s.xy[k][j-1][i] + H_s.xz[k][j-1][i] - H_s.xy[k][j][i] - H_s.xz[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                    Enp1 += Cc * E_nm1.zy[k][j][i] -
                            1. / 2. * Cb * params.delta.dy *
                                    ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dy * J_c.zy[k][j][i];

                if ((is_disp || params.is_disp_ml) && gamma_l) {
                    Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                            kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                    Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.zy[k][j][i];
                    E_nm1.zy[k][j][i] = E_s.zy[k][j][i];
                    J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                    J_s.zy[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + E_s.zy[k][j][i]); }

                eh_vec[n][j][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;
                ca_vec[n][j - 1] = Ca;
                cb_vec[n][j - 1] = Cb;
                }
                if (J_tot > 1) {
                j = 0;
                eh_vec[n][j][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;
                first_derivative(eh_vec[n], eh_vec[n], dk_e_y, N_e_y, E_s.zy.plan_f[n],
                                    E_s.zy.plan_b[n]);
                }
                for (j = 1; j < J_tot; j++) {
                E_s.zy[k][j][i] = ca_vec[n][j - 1] * E_s.zy[k][j][i] -
                                    cb_vec[n][j - 1] * eh_vec[n][j][0] / ((double) N_e_y);
                //E_s.zy[k][j][i] = Enp1;
                }
            }
            //PSTD, E_s.zy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
        }  //(params.dimension==THREE || params.dimension==TE)
        else {
    #pragma omp for
        for (k = 0; k <= K_tot; k++)
            for (j = 1; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
                rho = 0.;
                k_loc = k;
                if (params.is_structure)
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                    k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                    k_loc = params.pml.Dzl + K - 1;
                    else
                    k_loc = params.pml.Dzl + 1;
                }
                if (!params.is_multilayer) array_ind = j;
                else
                array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (!materials[k][j][i]) {
                Ca = C.a.y[array_ind];
                Cb = C.b.y[array_ind];
                if (params.is_disp_ml) Cc = C.c.y[array_ind];
                else
                    Cc = 0.;
                if (is_conductive) rho = rho_cond.y[array_ind];
                } else {
                rho = 0.;
                Ca = Cmaterial.a.y[materials[k][j][i] - 1];
                Cb = Cmaterial.b.y[materials[k][j][i] - 1];
                Cc = Cmaterial.c.y[materials[k][j][i] - 1];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || params.is_disp_ml) {
                kappa_l = matched_layer.kappa.y[array_ind];
                sigma_l = matched_layer.sigma.y[array_ind];
                alpha_l = matched_layer.alpha[k_loc];
                beta_l = matched_layer.beta[k_loc];
                gamma_l = matched_layer.gamma[k_loc];

                if (!materials[k][j][i]) {
                    alpha_l = 0.;
                    beta_l = 0.;
                    gamma_l = 0.;
                } else {
                    alpha_l = matched_layer.alpha[k_loc];
                    beta_l = matched_layer.beta[k_loc];
                    gamma_l = matched_layer.gamma[k_loc];
                }
                }


                Enp1 = Ca * E_s.zy[k][j][i] + Cb * (H_s.xy[k][j - 1][i] + H_s.xz[k][j - 1][i] -
                                                    H_s.xy[k][j][i] - H_s.xz[k][j][i]);
                if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zy[k][j][i] -
                        1. / 2. * Cb * params.delta.dy *
                                ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * params.delta.dy * J_c.zy[k][j][i];

                if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                        kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                Jnp1 += sigma_l / EPSILON0 * gamma_l * E_s.zy[k][j][i];
                E_nm1.zy[k][j][i] = E_s.zy[k][j][i];
                J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                J_s.zy[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + E_s.zy[k][j][i]); }

                E_s.zy[k][j][i] = Enp1;
            }
        }
    }//end of parallel section
    //fprintf(stderr,"Pos 09:\n");
    if (TIME_EXEC) {
        click_timer(IterationTimers::INTERNAL);
        ;
    }
    /********************/

    //update terms for self consistency across scattered/total interface - E updates##
    if (params.source_mode == SourceMode::steadystate) {//steadystate
        complex<double> commonPhase =
                exp(-IMAGINARY_UNIT * fmod(params.omega_an * time_H, 2. * DCPI));
        double commonAmplitude = linear_ramp(time_H, 1. / (params.omega_an / (2 * DCPI)));
        for (k = (K0.index); k <= (K1.index); k++)
        for (j = (J0.index); j <= (J1.index); j++) {
            if (I0.apply) {//Perform across I0

            if (!params.is_multilayer) array_ind = I0.index;
            else
                array_ind = (I_tot + 1) * k + I0.index;

            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
                E_s.zx[k][j][I0.index] =
                        E_s.zx[k][j][I0.index] -
                        C.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][2]));
                if (is_conductive)
                J_c.zx[k][j][I0.index] +=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
                if (params.is_disp_ml)
                J_s.zx[k][j][I0.index] +=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
            }
            if (j < (J1.index)) {
                E_s.yx[k][j][I0.index] =
                        E_s.yx[k][j][I0.index] +
                        C.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][3]));
                if (is_conductive)
                J_c.yx[k][j][I0.index] -=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
                if (params.is_disp_ml)
                J_s.yx[k][j][I0.index] -=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
            }
            }
            if (I1.apply) {//Perform across I1

            if (!params.is_multilayer) array_ind = I1.index;
            else
                array_ind = (I_tot + 1) * k + I1.index;

            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
                E_s.zx[k][j][I1.index] =
                        E_s.zx[k][j][I1.index] +
                        C.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][6]));
                if (is_conductive)
                J_c.zx[k][j][I1.index] -=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
                if (params.is_disp_ml)
                J_s.zx[k][j][I1.index] -=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
            }
            if (j < (J1.index)) {
                E_s.yx[k][j][I1.index] =
                        E_s.yx[k][j][I1.index] -
                        C.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][7]));
                if (is_conductive)
                J_c.yx[k][j][I1.index] +=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
                if (params.is_disp_ml)
                J_s.yx[k][j][I1.index] +=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                                IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
            }
            }
        }

        for (k = (K0.index); k <= (K1.index); k++)
        for (i = (I0.index); i <= (I1.index); i++) {
            if (J0.apply) {//Perform across J0
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC) {

                if (!params.is_multilayer) array_ind = J0.index;
                else
                array_ind = (J_tot + 1) * k + J0.index;

                E_s.zy[k][(J0.index)][i] =
                        E_s.zy[k][(J0.index)][i] +
                        C.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
                if (is_conductive)
                J_c.zy[k][(J0.index)][i] -=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
                if (params.is_disp_ml)
                J_s.zy[k][(J0.index)][i] -=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
            }
            if (i < (I1.index)) {
                E_s.xy[k][(J0.index)][i] =
                        E_s.xy[k][(J0.index)][i] -
                        C.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
                if (is_conductive)
                J_c.xy[k][(J0.index)][i] +=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
                if (params.is_disp_ml)
                J_s.xy[k][(J0.index)][i] +=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
            }
            }
            if (J1.apply) {//Perform across J1

            if (!params.is_multilayer) array_ind = J1.index;
            else
                array_ind = (J_tot + 1) * k + J1.index;

            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
                E_s.zy[k][(J1.index)][i] =
                        E_s.zy[k][(J1.index)][i] -
                        C.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
                if (is_conductive)
                J_c.zy[k][(J1.index)][i] +=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
                if (params.is_disp_ml)
                J_s.zy[k][(J1.index)][i] -=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
            }
            if (i < (I1.index)) {
                E_s.xy[k][(J1.index)][i] =
                        E_s.xy[k][(J1.index)][i] +
                        C.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
                if (is_conductive)
                J_c.xy[k][(J1.index)][i] -=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
                if (params.is_disp_ml)
                J_s.xy[k][(J1.index)][i] +=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                                (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                                IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
            }
            }
        }

        for (j = (J0.index); j <= (J1.index); j++)
        for (i = (I0.index); i <= (I1.index); i++) {
            if (K0.apply) {//Perform across K0
            if (j < (J1.index)) {
                E_s.yz[(K0.index)][j][i] =
                        E_s.yz[(K0.index)][j][i] -
                        C.b.z[K0.index] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
                if (is_conductive)
                J_c.yz[(K0.index)][j][i] +=
                        rho_cond.z[(K0.index)] * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
                if (params.is_disp_ml)
                J_s.yz[(K0.index)][j][i] -=
                        matched_layer.kappa.z[(K0.index)] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
            }
            if (i < (I1.index)) {
                E_s.xz[(K0.index)][j][i] =
                        E_s.xz[(K0.index)][j][i] +
                        C.b.z[K0.index] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
                if (is_conductive)
                J_c.xz[(K0.index)][j][i] -=
                        rho_cond.z[(K0.index)] * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
                if (params.is_disp_ml)
                J_s.xz[(K0.index)][j][i] +=
                        matched_layer.kappa.z[(K0.index)] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
            }
            }
            if (K1.apply) {//Perform across K1
            if (j < (J1.index)) {
                E_s.yz[(K1.index)][j][i] =
                        E_s.yz[(K1.index)][j][i] +
                        C.b.z[K1.index] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
                if (is_conductive)
                J_c.yz[(K1.index)][j][i] -=
                        rho_cond.z[(K1.index)] * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
                if (params.is_disp_ml)
                J_s.yz[(K1.index)][j][i] +=
                        matched_layer.kappa.z[(K1.index)] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
            }
            if (i < (I1.index)) {
                E_s.xz[(K1.index)][j][i] =
                        E_s.xz[(K1.index)][j][i] -
                        C.b.z[K1.index] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
                if (is_conductive)
                J_c.xz[(K1.index)][j][i] +=
                        rho_cond.z[(K1.index)] * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
                if (params.is_disp_ml)
                J_s.xz[(K1.index)][j][i] -=
                        matched_layer.kappa.z[(K1.index)] * matched_layer.gamma[k] /
                        (2. * params.dt) * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                                (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
            }
            }
        }
        H.ft = real(commonAmplitude * commonPhase);
    } else if (params.source_mode == SourceMode::pulsed) {//pulsed

        if (J_tot == 0) {
        j = 0;
        for (i = 0; i < (I_tot + 1); i++) {
            E_s.yz[K0.index][j][i] =
                    E_s.yz[K0.index][j][i] -
                    C.b.z[K0.index] *
                            real((Ksource.real[0][i - (I0.index)][2] +
                                IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][2]) *
                                (-1.0 * IMAGINARY_UNIT) *
                                exp(-IMAGINARY_UNIT *
                                    fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                            exp(-1.0 * DCPI *
                                pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) /
                                            (params.hwhm),
                                    2));
            //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (is_conductive)
            J_c.yz[K0.index][j][i] +=
                    rho_cond.z[K0.index] * C.b.z[K0.index] *
                    real((Ksource.real[0][i - (I0.index)][2] +
                            IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][2]) *
                            (-1.0 * IMAGINARY_UNIT) *
                            exp(-IMAGINARY_UNIT *
                                fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm),
                            2));
            //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (params.is_disp_ml) {
            J_s.yz[K0.index][j][i] -=
                    matched_layer.kappa.z[K0.index] * matched_layer.gamma[K0.index] /
                    (2. * params.dt) * C.b.z[K0.index] *
                    real((Ksource.real[0][i - (I0.index)][2] +
                            IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][2]) *
                            (-1.0 * IMAGINARY_UNIT) *
                            exp(-IMAGINARY_UNIT *
                                fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm),
                            2));
            //J_s.yz[(int)K0[0]][j][i] -= matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            }
        }
        } else
        for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
            /*
        if(i==41 & j==41)
        fprintf(stderr,"C.b.z = %.10e, Re(K) = %.10e, Im(K) = %.10e, time_H= %.10e, params.to_l=%.10e, params.delta.dz/LIGHT_V/2=%.10e, hwhm = %.10e, dE=%.10e\n",C.b.z[(int)K0[0]],Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2],Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2],time_H,params.to_l,params.delta.dz/LIGHT_V/2,params.hwhm,C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l + params.delta.dz/LIGHT_V/2.)/(params.hwhm),2)));
        */
            E_s.yz[K0.index][j][i] =
                    E_s.yz[K0.index][j][i] -
                    C.b.z[K0.index] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                                    (-1.0 * IMAGINARY_UNIT) *
                                    exp(-IMAGINARY_UNIT *
                                        fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                            exp(-1.0 * DCPI *
                                pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) /
                                            (params.hwhm),
                                    2));
            //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (is_conductive)
                J_c.yz[K0.index][j][i] +=
                        rho_cond.z[K0.index] * C.b.z[K0.index] *
                        real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                            IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                            (-1.0 * IMAGINARY_UNIT) *
                            exp(-IMAGINARY_UNIT *
                                fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                        exp(-1.0 * DCPI *
                            pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) /
                                        (params.hwhm),
                                2));
            //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (params.is_disp_ml) {
                J_s.yz[K0.index][j][i] -=
                        matched_layer.kappa.z[K0.index] * matched_layer.gamma[K0.index] /
                        (2. * params.dt) * C.b.z[K0.index] *
                        real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                            IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                            (-1.0 * IMAGINARY_UNIT) *
                            exp(-IMAGINARY_UNIT *
                                fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                        exp(-1.0 * DCPI *
                            pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) /
                                        (params.hwhm),
                                2));
                //J_s.yz[(int)K0[0]][j][i] -= matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            }
            }
        for (j = 0; j < (J_tot + 1); j++)
        for (i = 0; i < I_tot; i++) {
            E_s.xz[K0.index][j][i] =
                    E_s.xz[K0.index][j][i] +
                    C.b.z[K0.index] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                                IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                                (-1.0 * IMAGINARY_UNIT) *
                                exp(-IMAGINARY_UNIT *
                                    fmod(params.omega_an * (time_H - params.to_l), 2 * DCPI))) *
                            exp(-1.0 * DCPI *
                                pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) /
                                            (params.hwhm),
                                    2));
            //E_s.xz[(int)K0[0]][j][i] = E_s.xz[(int)K0[0]][j][i] + C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
            if (is_conductive)
            J_c.xz[K0.index][j][i] -=
                    rho_cond.z[K0.index] * C.b.z[K0.index] *
                    real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                            IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                            (-1.0 * IMAGINARY_UNIT) *
                            exp(-IMAGINARY_UNIT *
                                fmod(params.omega_an * (time_H - params.to_l), 2 * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm),
                            2));
            //J_c.xz[(int)K0[0]][j][i] -= rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
            if (params.is_disp_ml)
            J_s.xz[K0.index][j][i] +=
                    matched_layer.kappa.z[K0.index] * matched_layer.gamma[K0.index] /
                    (2. * params.dt) * C.b.z[K0.index] *
                    real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                            IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                            (-1.0 * IMAGINARY_UNIT) *
                            exp(-IMAGINARY_UNIT *
                                fmod(params.omega_an * (time_H - params.to_l), 2 * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm),
                            2));
            //J_s.xz[(int)K0[0]][j][i] += matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
        }
        //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
        H.ft = real((-1.0 * IMAGINARY_UNIT) *
                    exp(-IMAGINARY_UNIT *
                        fmod(params.omega_an * (time_H - params.to_l), 2. * DCPI))) *
                exp(-1.0 * DCPI *
                    pow((time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
        //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
    }
    //fprintf(stderr,"Pos 10:\n");

    //end of source terms
    if (TIME_EXEC) {
        click_timer(IterationTimers::INTERNAL);
        ;
    }

    /********************/
    //begin parallel
    #pragma omp parallel default(shared) private(i, j, k, n, k_loc,                                    \
                                                array_ind)//,ca_vec,cb_vec,cc_vec,eh_vec)
    {
        n = omp_get_thread_num();

        if (params.dimension == THREE || params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, H_s.xz
    #pragma omp for
            //H_s.xz updates
            for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_bound; j++)
                for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }

                if (!materials[k][j][i])
                    H_s.xz[k][j][i] = D.a.z[k_loc] * H_s.xz[k][j][i] +
                                    D.b.z[k_loc] * (E_s.yx[k + 1][j][i] + E_s.yz[k + 1][j][i] -
                                                    E_s.yx[k][j][i] - E_s.yz[k][j][i]);
                else
                    H_s.xz[k][j][i] = Dmaterial.a.z[materials[k][j][i] - 1] * H_s.xz[k][j][i] +
                                    Dmaterial.b.z[materials[k][j][i] - 1] *
                                            (E_s.yx[k + 1][j][i] + E_s.yz[k + 1][j][i] -
                                                E_s.yx[k][j][i] - E_s.yz[k][j][i]);
                }
            //FDTD, H_s.xz
        } else {
    #pragma omp for
            //H_s.xz updates
            for (j = 0; j < J_tot_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
                for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }

                if (!materials[k][j][i]) {
                    ca_vec[n][k] = D.a.z[k_loc];
                    cb_vec[n][k] = D.b.z[k_loc];
                    //H_s.xz[k][j][i] = D.a.z[k_loc]*H_s.xz[k][j][i]+D.b.z[k_loc]*(E_s.yx[k+1][j][i] + E_s.yz[k+1][j][i] - E_s.yx[k][j][i] - E_s.yz[k][j][i]);
                } else {
                    ca_vec[n][k] = Dmaterial.a.z[materials[k][j][i] - 1];
                    cb_vec[n][k] = Dmaterial.b.z[materials[k][j][i] - 1];
                    //H_s.xz[k][j][i] = Dmaterial.Da.z[materials[k][j][i]-1]*H_s.xz[k][j][i]+Dmaterial.Db.z[materials[k][j][i]-1]*(E_s.yx[k+1][j][i] + E_s.yz[k+1][j][i] - E_s.yx[k][j][i] - E_s.yz[k][j][i]);
                }

                eh_vec[n][k][0] = E_s.yx[k][j][i] + E_s.yz[k][j][i];
                eh_vec[n][k][1] = 0.;
                }
                k = K_tot;
                eh_vec[n][k][0] = E_s.yx[k][j][i] + E_s.yz[k][j][i];
                eh_vec[n][k][1] = 0.;

                first_derivative(eh_vec[n], eh_vec[n], dk_h_z, N_h_z, H_s.xz.plan_f[n],
                                H_s.xz.plan_b[n]);

                for (k = 0; k < K_tot; k++) {
                H_s.xz[k][j][i] = ca_vec[n][k] * H_s.xz[k][j][i] +
                                    cb_vec[n][k] * eh_vec[n][k][0] / ((double) N_h_z);
                }
            }

            //PSTD, H_s.xz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)

        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, H_s.xy
    #pragma omp for
            //H_s.xy updates
            for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot; j++)
                for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;
                if (!materials[k][j][i])
                    H_s.xy[k][j][i] = D.a.y[array_ind] * H_s.xy[k][j][i] +
                                    D.b.y[array_ind] * (E_s.zy[k][j][i] + E_s.zx[k][j][i] -
                                                        E_s.zy[k][j + 1][i] - E_s.zx[k][j + 1][i]);
                else
                    H_s.xy[k][j][i] = Dmaterial.a.y[materials[k][j][i] - 1] * H_s.xy[k][j][i] +
                                    Dmaterial.b.y[materials[k][j][i] - 1] *
                                            (E_s.zy[k][j][i] + E_s.zx[k][j][i] -
                                                E_s.zy[k][j + 1][i] - E_s.zx[k][j + 1][i]);
                }
            //FDTD, H_s.xy
        } else {
    #pragma omp for
            //H_s.xy updates
            for (k = 0; k < K_tot; k++)
            for (i = 0; i < (I_tot + 1); i++) {
                for (j = 0; j < J_tot; j++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;
                if (!materials[k][j][i]) {
                    ca_vec[n][j] = D.a.y[array_ind];
                    cb_vec[n][j] = D.b.y[array_ind];
                    //		H_s.xy[k][j][i] = D.a.y[array_ind]*H_s.xy[k][j][i]+D.b.y[array_ind]*(E_s.zy[k][j][i] + E_s.zx[k][j][i] - E_s.zy[k][j+1][i] - E_s.zx[k][j+1][i]);
                } else {
                    ca_vec[n][j] = Dmaterial.a.y[materials[k][j][i] - 1];
                    cb_vec[n][j] = Dmaterial.b.y[materials[k][j][i] - 1];
                    //		H_s.xy[k][j][i] = Dmaterial.Da.y[materials[k][j][i]-1]*H_s.xy[k][j][i]+Dmaterial.Db.y[materials[k][j][i]-1]*(E_s.zy[k][j][i] + E_s.zx[k][j][i] - E_s.zy[k][j+1][i] - E_s.zx[k][j+1][i]);
                }

                eh_vec[n][j][0] = E_s.zy[k][j][i] + E_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;
                }
                j = J_tot;
                eh_vec[n][j][0] = E_s.zy[k][j][i] + E_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;

                first_derivative(eh_vec[n], eh_vec[n], dk_h_y, N_h_y, H_s.xy.plan_f[n],
                                H_s.xy.plan_b[n]);

                for (j = 0; j < J_tot; j++) {
                H_s.xy[k][j][i] = ca_vec[n][j] * H_s.xy[k][j][i] -
                                    cb_vec[n][j] * eh_vec[n][j][0] / ((double) N_h_y);
                }

                /*
    if( i==12 && k==24){
    fprintf(stdout,"tind: %d\n",tind);
    fprintf(stdout,"Da: ");
    for(j=0;j<J_tot;j++)
    fprintf(stdout,"%e ",ca_vec[n][j]);
    fprintf(stdout,"\nDb: ");
    for(j=0;j<J_tot;j++)
    fprintf(stdout,"%e ",cb_vec[n][j]);

    fprintf(stdout,"\neh_vec: ");
    for(j=0;j<J_tot;j++)
    fprintf(stdout,"%e ",eh_vec[n][j][0]/((double) N_e_y));
    fprintf(stdout,"\n");
    }
        */
            }
            //PSTD, H_s.xy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)

        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, H_s.yx
    #pragma omp for
            //H_s.yx updates
            for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++)
                for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;
                if (!materials[k][j][i])
                    H_s.yx[k][j][i] = D.a.x[array_ind] * H_s.yx[k][j][i] +
                                    D.b.x[array_ind] * (E_s.zx[k][j][i + 1] + E_s.zy[k][j][i + 1] -
                                                        E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                else {
                    H_s.yx[k][j][i] = Dmaterial.a.x[materials[k][j][i] - 1] * H_s.yx[k][j][i] +
                                    Dmaterial.b.x[materials[k][j][i] - 1] *
                                            (E_s.zx[k][j][i + 1] + E_s.zy[k][j][i + 1] -
                                                E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                }
                }
            //FDTD, H_s.yx
        } else {
    #pragma omp for
            //H_s.yx updates
            for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++) {
                for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;
                if (!materials[k][j][i]) {
                    ca_vec[n][i] = D.a.x[array_ind];
                    cb_vec[n][i] = D.b.x[array_ind];
                    //		H_s.yx[k][j][i] = D.a.x[array_ind]*H_s.yx[k][j][i]+D.b.x[array_ind]*(E_s.zx[k][j][i+1] + E_s.zy[k][j][i+1] - E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                } else {
                    ca_vec[n][i] = Dmaterial.a.x[materials[k][j][i] - 1];
                    cb_vec[n][i] = Dmaterial.b.x[materials[k][j][i] - 1];
                    //	H_s.yx[k][j][i] = Dmaterial.Da.x[materials[k][j][i]-1]*H_s.yx[k][j][i]+Dmaterial.Db.x[materials[k][j][i]-1]*(E_s.zx[k][j][i+1] + E_s.zy[k][j][i+1] - E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                }

                eh_vec[n][i][0] = E_s.zx[k][j][i] + E_s.zy[k][j][i];
                eh_vec[n][i][1] = 0.;
                }
                i = I_tot;
                eh_vec[n][i][0] = E_s.zx[k][j][i] + E_s.zy[k][j][i];
                eh_vec[n][i][1] = 0.;

                first_derivative(eh_vec[n], eh_vec[n], dk_h_x, N_h_x, H_s.yx.plan_f[n],
                                H_s.yx.plan_b[n]);

                for (i = 0; i < I_tot; i++) {
                H_s.yx[k][j][i] = ca_vec[n][i] * H_s.yx[k][j][i] +
                                    cb_vec[n][i] * eh_vec[n][i][0] / ((double) N_h_x);
                }
            }
            //PSTD, H_s.yx
        }

        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, H_s.yz
    #pragma omp for
            //H_s.yz updates
            for (k = 0; k < K_tot; k++) {
            for (j = 0; j < J_tot_p1_bound; j++)
                for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!materials[k][j][i]) {
                    /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                    H_s.yz[k][j][i] = D.a.z[k_loc] * H_s.yz[k][j][i] +
                                    D.b.z[k_loc] * (E_s.xy[k][j][i] + E_s.xz[k][j][i] -
                                                    E_s.xy[k + 1][j][i] - E_s.xz[k + 1][j][i]);
                } else {
                    /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                    H_s.yz[k][j][i] = Dmaterial.a.z[materials[k][j][i] - 1] * H_s.yz[k][j][i] +
                                    Dmaterial.b.z[materials[k][j][i] - 1] *
                                            (E_s.xy[k][j][i] + E_s.xz[k][j][i] -
                                                E_s.xy[k + 1][j][i] - E_s.xz[k + 1][j][i]);
                }
                }
            }
            //FDTD, H_s.yz
        } else {
            //#pragma omp for
            //H_s.yz updates
            for (j = 0; j < J_tot_p1_bound; j++)
    #pragma omp for
            for (i = 0; i < I_tot; i++) {
                for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!materials[k][j][i]) {
                    ca_vec[n][k] = D.a.z[k_loc];
                    cb_vec[n][k] = D.b.z[k_loc];
                    /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                    //H_s.yz[k][j][i] = D.a.z[k_loc]*H_s.yz[k][j][i]+D.b.z[k_loc]*(E_s.xy[k][j][i] + E_s.xz[k][j][i] - E_s.xy[k+1][j][i] - E_s.xz[k+1][j][i]);
                } else {
                    ca_vec[n][k] = Dmaterial.a.z[materials[k][j][i] - 1];
                    cb_vec[n][k] = Dmaterial.b.z[materials[k][j][i] - 1];
                    /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                    //H_s.yz[k][j][i] = Dmaterial.Da.z[materials[k][j][i]-1]*H_s.yz[k][j][i]+Dmaterial.Db.z[materials[k][j][i]-1]*(E_s.xy[k][j][i] + E_s.xz[k][j][i] - E_s.xy[k+1][j][i] - E_s.xz[k+1][j][i]);
                }

                eh_vec[n][k][0] = E_s.xy[k][j][i] + E_s.xz[k][j][i];
                eh_vec[n][k][1] = 0.;
                }
                k = K_tot;
                eh_vec[n][k][0] = E_s.xy[k][j][i] + E_s.xz[k][j][i];
                eh_vec[n][k][1] = 0.;

                /*
    if( i==12 & j==12 ){
    for(k=0;k<K_tot;k++)
    fprintf(stdout,"%.10e ",eh_vec[n][k][0]);
    fprintf(stdout,"\n");
    }
        */

                first_derivative(eh_vec[n], eh_vec[n], dk_h_z, N_h_z, H_s.yz.plan_f[n],
                                H_s.yz.plan_b[n]);

                for (k = 0; k < K_tot; k++) {
                H_s.yz[k][j][i] = ca_vec[n][k] * H_s.yz[k][j][i] -
                                    cb_vec[n][k] * eh_vec[n][k][0] / ((double) N_h_z);
                }
            }
            //PSTD, H_s.yz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
        }  //(params.dimension==THREE || params.dimension==TE)
        else {

    #pragma omp for
        for (k = 0; k <= K_tot; k++)
            for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++)
                if (!materials[k][j][i]) H_s.xz[k][j][i] = 0.;
                else
                H_s.xz[k][j][i] = 0.;

    #pragma omp for
        //H_s.xy update
        for (k = 0; k <= K_tot; k++)
            for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (params.is_structure)
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                    k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                    k_loc = params.pml.Dzl + K - 1;
                    else
                    k_loc = params.pml.Dzl + 1;
                }
                if (!params.is_multilayer) array_ind = j;
                else
                array_ind = (J_tot + 1) * k_loc + j;
                if (!materials[k][j][i])
                H_s.xy[k][j][i] = D.a.y[array_ind] * H_s.xy[k][j][i] +
                                    D.b.y[array_ind] * (E_s.zy[k][j][i] + E_s.zx[k][j][i] -
                                                        E_s.zy[k][j + 1][i] - E_s.zx[k][j + 1][i]);
                else
                H_s.xy[k][j][i] = Dmaterial.a.y[materials[k][j][i] - 1] * H_s.xy[k][j][i] +
                                    Dmaterial.b.y[materials[k][j][i] - 1] *
                                            (E_s.zy[k][j][i] + E_s.zx[k][j][i] - E_s.zy[k][j + 1][i] -
                                            E_s.zx[k][j + 1][i]);
            }

    #pragma omp for
        //H_s.yx update
        for (k = 0; k <= K_tot; k++)
            for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                    k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                    k_loc = params.pml.Dzl + K - 1;
                    else
                    k_loc = params.pml.Dzl + 1;
                }
                if (!params.is_multilayer) array_ind = i;
                else
                array_ind = (I_tot + 1) * k_loc + i;
                if (!materials[k][j][i])
                H_s.yx[k][j][i] = D.a.x[array_ind] * H_s.yx[k][j][i] +
                                    D.b.x[array_ind] * (E_s.zx[k][j][i + 1] + E_s.zy[k][j][i + 1] -
                                                        E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                else
                H_s.yx[k][j][i] = Dmaterial.a.x[materials[k][j][i] - 1] * H_s.yx[k][j][i] +
                                    Dmaterial.b.x[materials[k][j][i] - 1] *
                                            (E_s.zx[k][j][i + 1] + E_s.zy[k][j][i + 1] -
                                            E_s.zx[k][j][i] - E_s.zy[k][j][i]);
            }

    #pragma omp for
        for (k = 0; k <= K_tot; k++) {
            for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++)
                if (!materials[k][j][i]) H_s.yz[k][j][i] = 0.;
                else
                H_s.yz[k][j][i] = 0.;
        }
        }

        if (params.dimension == THREE || params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, H_s.zy
    #pragma omp for
            //H_s.zy update
            for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot; j++)
                for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;
                if (!materials[k][j][i])
                    H_s.zy[k][j][i] = D.a.y[array_ind] * H_s.zy[k][j][i] +
                                    D.b.y[array_ind] * (E_s.xy[k][j + 1][i] + E_s.xz[k][j + 1][i] -
                                                        E_s.xy[k][j][i] - E_s.xz[k][j][i]);
                else
                    H_s.zy[k][j][i] = Dmaterial.a.y[materials[k][j][i] - 1] * H_s.zy[k][j][i] +
                                    Dmaterial.b.y[materials[k][j][i] - 1] *
                                            (E_s.xy[k][j + 1][i] + E_s.xz[k][j + 1][i] -
                                                E_s.xy[k][j][i] - E_s.xz[k][j][i]);
                }
            //FDTD, H_s.zy
        } else {
    #pragma omp for
            //H_s.zy update
            for (k = 0; k < (K_tot + 1); k++)
            for (i = 0; i < I_tot; i++) {
                for (j = 0; j < J_tot; j++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = j;
                else
                    array_ind = (J_tot + 1) * k_loc + j;
                if (!materials[k][j][i]) {
                    ca_vec[n][j] = D.a.y[array_ind];
                    cb_vec[n][j] = D.b.y[array_ind];
                    //	      H_s.zy[k][j][i] = D.a.y[array_ind]*H_s.zy[k][j][i]+D.b.y[array_ind]*(E_s.xy[k][j+1][i] + E_s.xz[k][j+1][i] - E_s.xy[k][j][i] - E_s.xz[k][j][i]);
                } else {
                    ca_vec[n][j] = Dmaterial.a.y[materials[k][j][i] - 1];
                    cb_vec[n][j] = Dmaterial.b.y[materials[k][j][i] - 1];
                    //	      H_s.zy[k][j][i] = Dmaterial.Da.y[materials[k][j][i]-1]*H_s.zy[k][j][i]+Dmaterial.Db.y[materials[k][j][i]-1]*(E_s.xy[k][j+1][i] + E_s.xz[k][j+1][i] - E_s.xy[k][j][i] - E_s.xz[k][j][i]);
                }

                eh_vec[n][j][0] = E_s.xy[k][j][i] + E_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;
                }
                j = J_tot;
                eh_vec[n][j][0] = E_s.xy[k][j][i] + E_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;

                first_derivative(eh_vec[n], eh_vec[n], dk_h_y, N_h_y, H_s.zy.plan_f[n],
                                H_s.zy.plan_b[n]);

                for (j = 0; j < J_tot; j++) {
                H_s.zy[k][j][i] = ca_vec[n][j] * H_s.zy[k][j][i] +
                                    cb_vec[n][j] * eh_vec[n][j][0] / ((double) N_h_y);
                }
            }
            //PSTD, H_s.zy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)


        if (solver_method == SolverMethod::FiniteDifference) {
    //FDTD, H_s.zx
    #pragma omp for
            //H_s.zx update
            for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot_bound; j++)
                for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;
                if (!materials[k][j][i])
                    H_s.zx[k][j][i] = D.a.x[array_ind] * H_s.zx[k][j][i] +
                                    D.b.x[array_ind] * (E_s.yx[k][j][i] + E_s.yz[k][j][i] -
                                                        E_s.yx[k][j][i + 1] - E_s.yz[k][j][i + 1]);
                else
                    H_s.zx[k][j][i] = Dmaterial.a.x[materials[k][j][i] - 1] * H_s.zx[k][j][i] +
                                    Dmaterial.b.x[materials[k][j][i] - 1] *
                                            (E_s.yx[k][j][i] + E_s.yz[k][j][i] -
                                                E_s.yx[k][j][i + 1] - E_s.yz[k][j][i + 1]);
                }
            //FDTD, H_s.zx
        } else {
    #pragma omp for
            //H_s.zx update
            for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot_bound; j++) {
                for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (params.is_structure)
                    if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                    if ((k - structure[i][1]) < (K + params.pml.Dzl) &&
                        (k - structure[i][1]) > params.pml.Dzl)
                        k_loc = k - structure[i][1];
                    else if ((k - structure[i][1]) >= (K + params.pml.Dzl))
                        k_loc = params.pml.Dzl + K - 1;
                    else
                        k_loc = params.pml.Dzl + 1;
                    }
                if (!params.is_multilayer) array_ind = i;
                else
                    array_ind = (I_tot + 1) * k_loc + i;
                if (!materials[k][j][i]) {
                    //		H_s.zx[k][j][i] = D.a.x[array_ind]*H_s.zx[k][j][i]+D.b.x[array_ind]*(E_s.yx[k][j][i] + E_s.yz[k][j][i] - E_s.yx[k][j][i+1] - E_s.yz[k][j][i+1]);
                    ca_vec[n][i] = D.a.x[array_ind];
                    cb_vec[n][i] = D.b.x[array_ind];
                } else {
                    //		H_s.zx[k][j][i] = Dmaterial.Da.x[materials[k][j][i]-1]*H_s.zx[k][j][i]+Dmaterial.Db.x[materials[k][j][i]-1]*(E_s.yx[k][j][i] + E_s.yz[k][j][i] - E_s.yx[k][j][i+1] - E_s.yz[k][j][i+1]);
                    ca_vec[n][i] = Dmaterial.a.x[materials[k][j][i] - 1];
                    cb_vec[n][i] = Dmaterial.b.x[materials[k][j][i] - 1];
                }

                eh_vec[n][i][0] = E_s.yx[k][j][i] + E_s.yz[k][j][i];
                eh_vec[n][i][1] = 0.;
                }
                i = I_tot;
                eh_vec[n][i][0] = E_s.yx[k][j][i] + E_s.yz[k][j][i];
                eh_vec[n][i][1] = 0.;


                first_derivative(eh_vec[n], eh_vec[n], dk_h_x, N_h_x, H_s.zx.plan_f[n],
                                H_s.zx.plan_b[n]);

                for (i = 0; i < I_tot; i++) {
                H_s.zx[k][j][i] = ca_vec[n][i] * H_s.zx[k][j][i] -
                                    cb_vec[n][i] * eh_vec[n][i][0] / ((double) N_h_x);
                }
            }
            //PSTD, H_s.zx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
        }  //(params.dimension==THREE || params.dimension==TE)
    }    //end parallel
    if (TIME_EXEC) {
        click_timer(IterationTimers::INTERNAL);
        ;
    }

    //fprintf(stderr,"Pos 11b:\n");
    //update terms for self consistency across scattered/total interface - E updates
    if (params.source_mode == SourceMode::steadystate) {//steadystate
        complex<double> commonPhase =
                exp(-IMAGINARY_UNIT * fmod(params.omega_an * time_E, 2. * DCPI));
        double commonAmplitude = linear_ramp(time_E, 1. / (params.omega_an / (2 * DCPI)));
        for (k = (K0.index); k <= (K1.index); k++)
        for (j = (J0.index); j <= (J1.index); j++) {
            if (I0.apply) {//Perform across I0

            if (!params.is_multilayer) array_ind = I0.index - 1;
            else
                array_ind = (I_tot + 1) * k + I0.index - 1;

            if (j < (J1.index))
                H_s.zx[k][j][(I0.index) - 1] =
                        H_s.zx[k][j][(I0.index) - 1] +
                        D.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][0] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][0]));
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
                H_s.yx[k][j][(I0.index) - 1] =
                        H_s.yx[k][j][(I0.index) - 1] -
                        D.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][1] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][1]));
            }
            if (I1.apply) {//Perform across I1

            if (!params.is_multilayer) array_ind = I1.index;
            else
                array_ind = (I_tot + 1) * k + I1.index;

            if (j < (J1.index))
                H_s.zx[k][j][(I1.index)] =
                        H_s.zx[k][j][(I1.index)] -
                        D.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][4] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][4]));
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
                H_s.yx[k][j][(I1.index)] =
                        H_s.yx[k][j][(I1.index)] +
                        D.b.x[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Isource.real[k - (K0.index)][j - (J0.index)][5] +
                                    IMAGINARY_UNIT *
                                            Isource.imag[k - (K0.index)][j - (J0.index)][5]));
            }
        }

        for (k = (K0.index); k <= (K1.index); k++)
        for (i = (I0.index); i <= (I1.index); i++) {
            if (J0.apply) {//Perform across J0

            if (!params.is_multilayer) array_ind = J0.index;
            else
                array_ind = (J_tot + 1) * k + J0.index;

            if (i < (I1.index))
                H_s.zy[k][(J0.index) - 1][i] =
                        H_s.zy[k][(J0.index) - 1][i] -
                        D.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][0] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][0]));

            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
                H_s.xy[k][(J0.index) - 1][i] =
                        H_s.xy[k][(J0.index) - 1][i] +
                        D.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][1] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][1]));
            }
            if (J1.apply) {//Perform across J1

            if (!params.is_multilayer) array_ind = J1.index;
            else
                array_ind = (J_tot + 1) * k + J1.index;

            if (i < (I1.index))
                H_s.zy[k][(J1.index)][i] =
                        H_s.zy[k][(J1.index)][i] +
                        D.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][4] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][4]));
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
                H_s.xy[k][(J1.index)][i] =
                        H_s.xy[k][(J1.index)][i] -
                        D.b.y[array_ind] *
                                real(commonAmplitude * commonPhase *
                                    (Jsource.real[k - (K0.index)][i - (I0.index)][5] +
                                    IMAGINARY_UNIT *
                                            Jsource.imag[k - (K0.index)][i - (I0.index)][5]));
            }
        }

        for (j = (J0.index); j <= (J1.index); j++)
        for (i = (I0.index); i <= (I1.index); i++) {
            if (K0.apply) {//Perform across K0
            if (i < (I1.index))
                H_s.yz[(K0.index) - 1][j][i] =
                        H_s.yz[(K0.index) - 1][j][i] +
                        D.b.z[(K0.index) - 1] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][0] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][0]));
            if (j < (J1.index))
                H_s.xz[(K0.index) - 1][j][i] =
                        H_s.xz[(K0.index) - 1][j][i] -
                        D.b.z[(K0.index) - 1] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][1] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][1]));
            }
            if (K1.apply) {//Perform across K1
            if (i < (I1.index))
                H_s.yz[(K1.index)][j][i] =
                        H_s.yz[(K1.index)][j][i] -
                        D.b.z[(K1.index)] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][4] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][4]));
            if (j < (J1.index))
                H_s.xz[(K1.index)][j][i] =
                        H_s.xz[(K1.index)][j][i] +
                        D.b.z[(K1.index)] *
                                real(commonAmplitude * commonPhase *
                                    (Ksource.real[j - (J0.index)][i - (I0.index)][5] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][5]));
            }
        }
        E.ft = real(commonAmplitude * commonPhase);
    } else if (params.source_mode == 1) {//pulsed
        //fprintf(stderr,"Pos 11c\n");
        if (J_tot == 0) {
        //fprintf(stderr,"Pos 11d\n");
        j = 0;
        for (i = 0; i < (I_tot + 1); i++) {
            H_s.xz[(K0.index) - 1][j][i] =
                    H_s.xz[(K0.index) - 1][j][i] -
                    D.b.z[(K0.index) - 1] *
                            real((Ksource.real[0][i - (I0.index)][1] +
                                IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][1]) *
                                (-1. * IMAGINARY_UNIT) *
                                exp(-IMAGINARY_UNIT *
                                    fmod(params.omega_an * (time_E - params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.eyi_present)
            H_s.xz[(K0.index) - 1][j][i] =
                    H_s.xz[(K0.index) - 1][j][i] - D.b.z[(K0.index) - 1] * Ei.y[tind][j][i];
        }
        //fprintf(stderr,"Pos 11e\n");
        for (i = 0; i < I_tot; i++) {
            H_s.yz[(K0.index) - 1][j][i] =
                    H_s.yz[(K0.index) - 1][j][i] +
                    D.b.z[(K0.index) - 1] *
                            real((Ksource.real[0][i - (I0.index)][0] +
                                IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][0]) *
                                (-1. * IMAGINARY_UNIT) *
                                exp(-IMAGINARY_UNIT *
                                    fmod(params.omega_an * (time_E - params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.exi_present)
            H_s.yz[(K0.index) - 1][j][i] =
                    H_s.yz[(K0.index) - 1][j][i] + D.b.z[(K0.index) - 1] * Ei.x[tind][j][i];
            //if(i==511)
            //  fprintf(stdout,"%e\n",D.b.z[((int)K0[0])-1]*exi[tind][j][i]);
        }
        //fprintf(stderr,"Pos 11f\n");
        } else {
        //fprintf(stderr,"Pos 11g\n");
        for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
            H_s.xz[(K0.index) - 1][j][i] =
                    H_s.xz[(K0.index) - 1][j][i] -
                    D.b.z[(K0.index) - 1] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][1] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][1]) *
                                    (-1. * IMAGINARY_UNIT) *
                                    exp(-IMAGINARY_UNIT *
                                        fmod(params.omega_an * (time_E - params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.eyi_present)
                H_s.xz[(K0.index) - 1][j][i] =
                        H_s.xz[(K0.index) - 1][j][i] - D.b.z[(K0.index) - 1] * Ei.y[tind][j][i];
            }
        //fprintf(stderr,"Pos 11h\n");
        for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++) {
            H_s.yz[(K0.index) - 1][j][i] =
                    H_s.yz[(K0.index) - 1][j][i] +
                    D.b.z[(K0.index) - 1] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][0] +
                                    IMAGINARY_UNIT *
                                            Ksource.imag[j - (J0.index)][i - (I0.index)][0]) *
                                    (-1. * IMAGINARY_UNIT) *
                                    exp(-IMAGINARY_UNIT *
                                        fmod(params.omega_an * (time_E - params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.exi_present)
                H_s.yz[(K0.index) - 1][j][i] =
                        H_s.yz[(K0.index) - 1][j][i] + D.b.z[(K0.index) - 1] * Ei.x[tind][j][i];
            }
        //fprintf(stderr,"Pos 11i\n");
        }
        E.ft = real((-1. * IMAGINARY_UNIT) *
                    exp(-IMAGINARY_UNIT * fmod(params.omega_an * (time_E - params.to_l), 2 * DCPI))) *
                exp(-1. * DCPI * pow((time_E - params.to_l) / (params.hwhm), 2.));
        //fprintf(stderr,"Pos 11j\n");
    }
    if (TIME_EXEC) {
        click_timer(IterationTimers::INTERNAL);
        ;
    }

    if (params.exphasorssurface || params.exphasorsvolume || params.exdetintegral ||
        vertex_phasors.there_are_vertices_to_extract_at()) {
        if (params.source_mode == SourceMode::steadystate) {
        /*
        Each time a new acquisition period of harmonic illumination begins, all complex amplitudes
        (volume, surface etc.) are set back to 0 since the discrete Fourier transforms used to acquire
        these complex amplitudes starts again. In particular, the returned complex amplitudes will have
        been acquired during a single acquisition period of harmonic illumination. Note that, as
        explained above, the acquisition period is actually three periods of the harmonic waves
        fundamental period. The complex amplitudes are reset to 0 using calls such as:

        initialiseDouble3DArray(ExR, dims[0], dims[1], dims[2]);

    However, the normalisation factors are reset to 0 here.
        */

        if ((tind % Nsteps) == 0) {
            E.angular_norm = 0.0;
            H.angular_norm = 0.0;

            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            E_norm[ifx] = 0.;
            H_norm[ifx] = 0.;
            }
        }

        /*In the calls below, the following two lines of code are equivalent up to numerical precision:

        E.add_to_angular_norm(fte, tind, Nsteps, params);
        E.add_to_angular_norm(fte, tind % Nsteps, Nsteps, params);

        To understand why, first consult the lines of code above:

        Nsteps_tmp = ceil(2.*DCPI/omega_an[0]/dt[0]*3);
        dt[0] = 2.*DCPI/omega_an[0]*3/Nsteps_tmp;
        Nsteps = (int)lround(Nsteps_tmp);

        Where dt and Nsteps are set. The reason for the factor of 3 is that we will perform complex
        amplitude extraction over 3 fundamental periods of the monochromatic source. We can then make
        the following statement:

        T/dt*3=1/(f*dt)*3=Nsteps

        where T and f (omega=2*pi*f) are the period and frequency of the monochromatic source, respectively.

        Then consider the argument of the exponentional function on phasor_norm, called by add_to_angular_norm, where tind=n is used:

        i*omega*((double) (n+1))*dt (where fmod(.,2*DCPI) is ignored since this will not affect the result)

        The argument of this function simplifies to:

        i*omega*(tind+1)*dt=i*2*pi*f*(tind+1)*dt=i*2*pi*(tind+1)*3/Nsteps (using f*dt=3/Nsteps)

        Then, without loss of generallity, let tind = p*Nsteps + q, substituting into the above

        i*2*pi*(tind+1)*3/Nsteps = i*2*pi*(p*Nsteps + q)*3/Nsteps = i*2*pi*3*p + i*2*pi*q*3/Nsteps

        In which case exp(i*2*pi*3*p + i*2*pi*q*3/Nsteps) = exp(i*2*pi*q*3/Nsteps)

        If instead we use tind % Nsteps=n, we see that n=q, leading to the same exponential function as
        above. So the two cases are equivalent.
        */

        E.add_to_angular_norm(tind, Nsteps, params);
        H.add_to_angular_norm(tind, Nsteps, params);

        extract_phasor_norms_at_all_frequencies(Nsteps);
        } else {
        if ((tind - params.start_tind) % params.Np == 0) {

            E.add_to_angular_norm(tind, params.Npe, params);
            H.add_to_angular_norm(tind, params.Npe, params);

            extract_phasor_norms_at_all_frequencies(params.Npe);
        }
        }
    }
    if (TIME_EXEC) { click_timer(IterationTimers::INTERNAL); }

    log_current_iteration_and_max_fields();
    //fprintf(stderr,"Post-iter 3\n");
    if ((params.source_mode == SourceMode::steadystate) && (tind == (params.Nt - 1)) &&
        (params.run_mode == RunMode::complete) && params.exphasorsvolume) {
        fprintf(stdout, "Iteration limit reached, setting output fields to last complete DFT\n");
        E.set_values_from(E_copy);
    }
    //fprintf(stderr,"Post-iter 4\n");
    fflush(stdout);
    //fprintf(stderr,"Post-iter 5\n");
    //fprintf(stderr,"%s %d %d\n",tdfdirstr, strcmp(tdfdirstr,""),are_equal(tdfdirstr,""));
    if (params.has_tdfdir && (tind % params.Np) == 0) {
        fprintf(stderr, "Saving field\n");
        ex_td_field_exporter.export_field(E_s, skip_tdf, tind);
    }
    //fprintf(stderr,"Post-iter 6\n");
    /*write out fdtdgrid to a file*/
    /*
        MATFile *toutfile;
        char toutputfilename[100];
        if(tind % params.Np == 0){
        //if(tind <= 1000){
        sprintf(toutputfilename,"tdata/fdtdgrid_%04d.mat",tind);
        toutfile = matOpen(toutputfilename, "w");
        matPutVariable(toutfile, "fdtdgrid", (mxArray *)in_matrices[0]);
        matClose(toutfile);
        }
    */
    /*write out fdtdgrid to a file*/

    }//end of main iteration loop

    if (TIME_MAIN_LOOP) {
        end_timer(IterationTimers::MAIN);
        spdlog::info("Time (seconds) ellapsed in main loop: {0:e}",
                    time_ellapsed_by(IterationTimers::MAIN));
    }
}
