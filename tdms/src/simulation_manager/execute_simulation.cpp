#include "simulation_manager.h"

#include <omp.h>

#include "numerical_derivative.h"

using namespace std;
using namespace tdms_math_constants;
using namespace tdms_phys_constants;

// whether or not to time subprocess execution within the main loop
#define TIME_EXEC false
// whether to time the main loop
#define TIME_MAIN_LOOP true
// convergence threshold; terminates the steady state iterations
#define TOL 1e-6

/*This function will take in the following arguments and perform the
 entire simulation

  fdtdgrid
  Cmaterial
  Dmaterial
  C
  D
  freespace
  disp_params
  delta
  interface
  Isource
  Jsource
  Ksource
  grid_labels
  omega_an
  to_l
  hwhm
  Dxl
  Dxu
  Dyl
  Dyu
  Dzl
  Dzu
  Nt
  dt
  start_tind
  sourcemode
  runmode
  exphasorsvolume
  exphasorssurface
  intphasorsurface
  phasorsurface
  phasorinc
  dimension
  conductive_aux
  dispersive_aux
  structure
  tdfield

  fdtdgrid - A structure with the following members, each of which is a 3 dimensional
  array:

  fdtdgrid.Exy       (double)
  fdtdgrid.Exz
  fdtdgrid.Eyx
  fdtdgrid.Eyz
  fdtdgrid.Ezx
  fdtdgrid.Ezy
  fdtdgrid.Hxy
  fdtdgrid.Hxz
  fdtdgrid.Hyx
  fdtdgrid.Hyz
  fdtdgrid.Hzx
  fdtdgrid.Hzy

  fdtdgrid.materials (uint8)

  Cmaterial - A structure with the following members, each of which is a 1 dimensional
  array indexed by indices od scattering materials:

  Cmaterial.Cax      (double)
  Cmaterial.Cay
  Cmaterial.Caz
  Cmaterial.Cbx
  Cmaterial.Cby
  Cmaterial.Cbz
  Cmaterial.Ccx
  Cmaterial.Ccy
  Cmaterial.Ccz

  Dmaterial - A structure with the following members, each of which is a 1 dimensional
  array indexed by indices od scattering materials:

  Dmaterial.Cax      (double)
  Dmaterial.Cay
  Dmaterial.Caz
  Dmaterial.Cbx
  Dmaterial.Cby
  Dmaterial.Cbz

  C - A structure with the same elements as Cmaterial whose elements relate to the background
  region. Elements are not restricted to 1 dimension.

  D - A structure with the same elements as Cmaterial whose elements relate to the background
  region. Elements are not restricted to 1 dimension.

  freespace - A structure with the following members each of which is a double

  freespace.Cbx      (double)
  freespace.Cby
  freespace.Cbz
  freespace.Dbx
  freespace.Dby
  freespace.Dbz

  disp_params - A structure with the following members:

  disp_params.alpha      (double)
  disp_params.beta
  disp_params.gamma

  delta - A structure with the following elements representing Yee cell dimension

  delta.x (double)
  delta.y
  delta.z

  interface - A structure with the following members each of which is a double

  interface.I0 - [I0 apply]
  interface.I1 - [I1 apply]
  interface.J0 - [J0 apply]
  interface.J1 - [J1 apply]
  interface.K0 - [K0 apply]
  interface.K1 - [K1 apply]

  where, for example, I0 is the position of the interface plane between
  total and scattered formulations. apply is a boolean which indicates
  whether or not to apply the boundary condition

  Isource - A 3d matrix of dimension 8x(J1-J0+1)x(K1-K0+1) which contains
  complex amplitude information for the interface equations at the I0 and I1
  planes. Should be Isource(1,:,:) - Ey, Isource(2,:,:) - Ez, Isource(3,:,:) - Hy,
  Isource(4,:,:) - Hz,Isource(5,:,:) - Ey, Isource(6,:,:) - Ez, Isource(7,:,:) - Hy,
  Isource(8,:,:) - Hz

  Jsource - A 3d matrix of dimension 8x(I1-I0+1)x(K1-K0+1) which contains
  complex amplitude information for the interface equations at the J0 and J1
  planes. Should be Jsource(1,:,:) - Ex, Jsource(2,:,:) - Ez, Jsource(3,:,:) - Hx,
  Jsource(4,:,:) - Hz,Jsource(5,:,:) - Ex, Jsource(6,:,:) - Ez, Jsource(7,:,:) - Hx,
  Jsource(8,:,:) - Hz

  Ksource - A 3d matrix of dimension 8x(I1-I0+1)x(J1-J0+1) which contains
  complex amplitude information for the interface equations at the K0 and K1
  planes. Should be Ksource(1,:,:) - Ex, Ksource(2,:,:) - Ey, Ksource(3,:,:) - Hx,
  Ksource(4,:,:) - Hy,Ksource(5,:,:) - Ex, Ksource(6,:,:) - Ey, Ksource(7,:,:) - Hx,
  Ksource(8,:,:) - Hy

  grid_labels - A structure with 3 elements, represents the axis labels:
  x_grid_label - cartesian labels of Yee cell origins (x-coordinates)
  y_grid_label - cartesian labels of Yee cell origins (y-coordinates)
  z_grid_label - cartesian labels of Yee cell origins (z-coordinates)

  omega_an - analytic angular frequency source

  to_l - time delay of pulse

  hwhm - hwhm of pulse

  Dxl - thickness of upper pml in the x direction
  Dxu - thickness of upper pml in the x direction
  Dyl - thickness of upper pml in the y direction
  Dyu - thickness of upper pml in the y direction
  Dzl - thickness of upper pml in the z direction
  Dzu - thickness of upper pml in the z direction

  lower_boundary_update - boolean for applying the lower update equations

  Nt - the number of iterations

  start_tind - the starting value of the tind index which controls the current time
  step of each fdtd iteration. The default for this would be 0;

  sourcemode - integer of value 0 for s steady state simulation and value 1 for a pulse simulation

  runmode - integer of value 0 for "complete" and 1 for "analyse"

  exphasorsvolume - Extract phasors in the FDTD volume if set to true

  exphasorssurface - Extract phasors about a specified cuboid surface if set to true

  intphasorssurface - Interpolate surface phasors onto a common point if true

  phasorsurface - A list of indices defining the cuboid to extract the phasors at

  phasorinc - An integer vector of three elements describing the factor by which to reduce the
  density of vertices in the enclosing observation mesh

  dimension - A string of value "3", "TE" or "TM"

  conductive_aux - auxiliary parameters required to model conductive multilayer

  dispersive_aux - auxiliary parameters required to model dispersive multilayer

  structure - 2 x (I_tot+1) integer array describing the grating structure, if one is present

  f_ex_vec - 1xN or Nx1 vector of frequencies at which to perform the extraction of complex amplitudes

  tdfield - structure containing elements exi and eyi which have dimension (I1-I0+1)x(J1-J0+1)xNt

  fieldsample.i - indices along the x-direction of locations at which to sample the field
  fieldsample.j - indices along the y-direction of locations at which to sample the field
  fieldsample.k - indices along the z-direction of locations at which to sample the field
  fieldsample.n - vector of the moments of the field to sample

  campssample.vertices - N x 3 matrix of indices where the complex amplitudes will be sampled
  campssample.components - numerical array of up to six elements which defines which field components
                           will be sampled, 1 means Ex, 2 Ey etc.
*/

void SimulationManager::execute_simulation() {
  spdlog::info("Using {} OMP threads", omp_get_max_threads());

  double rho;
  double alpha_l, beta_l, gamma_l;
  double kappa_l, sigma_l;
  double t0;

  double Ca, Cb, Cc;//used by interpolation scheme
  //the C and D vars for free space and pml
  double Enp1, Jnp1;

  double maxfield = 0.;

  double phaseTermE;
  complex<double> cphaseTermE;
  double lambda_an_t;

  int i, j, k, n, k_loc;
  int dft_counter = 0;

  complex<double> Idxt, Idyt, kprop;

  // get the number of Yee cells in each axial direction
  IJKDims IJK_tot = n_Yee_cells();
  int I_tot = IJK_tot.I_tot(), J_tot = IJK_tot.J_tot(), K_tot = IJK_tot.K_tot();

  // setup the variables that are to be used in the loop
  LoopVariables loop_variables(inputs, outputs.get_E_dimensions());

  /*Start of FDTD iteration*/
  //open a file for logging the times
  /*The times of the E and H fields at the point where update equations are applied.
    time_H is actually the time of the H field when the E field consistency update is
    applied and vice versa. time_E > time_H below since after the E field consistency
    update the E field will have advanced one time step.

    The interpretation of time is slightly complicated in the following. In what follows
    I write (tind*dt,(tind+1/2)*dt) to mean the time at which we currently know the
    electric (tind*dt) and magnetic ( (tind+1/2)*dt ) fields.

    Times before                Operation         Times after
    (tind*dt,(tind+1/2)*dt)     Extract phasors   (tind*dt,(tind+1/2)*dt)
    (tind*dt,(tind+1/2)*dt)     E field update    ( (tind+1)*dt,(tind+1/2)*dt)
    ((tind+1)*dt,(tind+1/2)*dt) H field update    ( (tind+1)*dt,(tind+3/2)*dt)
    ((tind+1)*dt,(tind+3/2)*dt) Normalisation extraction

    We note that the extractPhasorENorm uses (tind+1)*dt and extractPhasorHNorm uses
    (tind+1/2)*dt to perform the update equation in the DFT. This seems incorrect
    at first but we note that they take the terms fte and fth as inputs respectively.
    When one notes that fte is calculated using time_E and fth using time_H we see
    that this indexing is correct, ie, time_E = (tind+1)*dt and time_H = (tind+1/2)*dt.
  */
  double time_E;//< The real-"time" that the E field is currently at
  double time_H;//< The real-"time" that the H field is currently at
  t0 = (double) time(NULL);

  spdlog::debug("Starting main loop");

  if (TIME_MAIN_LOOP) { timers.start_timer(TimersTrackingLoop::MAIN); }

  for (unsigned int tind = inputs.params.start_tind; tind < inputs.params.Nt; tind++) {
    //fprintf(stderr,"Pos 00:\n");
    time_E = ((double) (tind + 1)) * inputs.params.dt;
    time_H = time_E - inputs.params.dt / 2.;
    //Extract phasors
    if ((dft_counter == inputs.Nsteps) && (inputs.params.run_mode == RunMode::complete) &&
        (inputs.params.source_mode == SourceMode::steadystate) && inputs.params.exphasorsvolume) {

      dft_counter = 0;

      double tol = outputs.E.normalised_difference(loop_variables.E_copy);
      if (tol < TOL) break;//required accuracy obtained

      spdlog::debug("Phasor convergence: {} (actual) > {} (required)", tol, TOL);
      loop_variables.E_copy.set_values_from(outputs.E);

      outputs.E.zero();
      outputs.H.zero();
      spdlog::debug("Zeroed the phasors");

      if (inputs.params.exphasorssurface) {
        outputs.surface_phasors.zero_surface_EH();
        spdlog::debug("Zeroed the surface components");
      }
    }

    if ((inputs.params.source_mode == SourceMode::steadystate) &&
        (inputs.params.run_mode == RunMode::complete) && inputs.params.exphasorsvolume) {

      outputs.E.set_phasors(inputs.E_s, dft_counter - 1, inputs.params.omega_an, inputs.params.dt,
                            inputs.Nsteps);
      outputs.H.set_phasors(inputs.H_s, dft_counter, inputs.params.omega_an, inputs.params.dt,
                            inputs.Nsteps);

      if (inputs.params.exphasorssurface) {
        for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
          outputs.surface_phasors.extractPhasorsSurface(
                  ifx, inputs.E_s, inputs.H_s, dft_counter, inputs.f_ex_vec[ifx] * 2 * DCPI,
                  inputs.Nsteps, inputs.params, inputs.params.intphasorssurface);
        }
        dft_counter++;
      }

    } else if ((inputs.params.source_mode == SourceMode::pulsed) &&
               (inputs.params.run_mode == RunMode::complete) && inputs.params.exphasorsvolume) {
      if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        outputs.E.set_phasors(inputs.E_s, tind - 1, inputs.params.omega_an, inputs.params.dt,
                              inputs.params.Npe);
        outputs.H.set_phasors(inputs.H_s, tind, inputs.params.omega_an, inputs.params.dt,
                              inputs.params.Npe);
      }
      if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
      //fprintf(stderr,"Pos 01b:\n");
    }

    /*extract fieldsample*/
    if (outputs.fieldsample.all_vectors_are_non_empty()) {
      outputs.fieldsample.extract(inputs.E_s, inputs.params.pml, inputs.params.Nt);
    }
    /*end extract fieldsample*/

    //fprintf(stderr,"Pos 02:\n");
    if (inputs.params.source_mode == SourceMode::pulsed &&
        inputs.params.run_mode == RunMode::complete && inputs.params.exphasorssurface) {
      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
          outputs.surface_phasors.extractPhasorsSurface(
                  ifx, inputs.E_s, inputs.H_s, tind, inputs.f_ex_vec[ifx] * 2 * DCPI,
                  inputs.params.Npe, inputs.params, inputs.params.intphasorssurface);
        }
      }
    }

    if (inputs.params.source_mode == SourceMode::pulsed &&
        inputs.params.run_mode == RunMode::complete &&
        (outputs.vertex_phasors.there_are_vertices_to_extract_at()) &&
        ((tind - inputs.params.start_tind) % inputs.params.Np == 0)) {
      //     fprintf(stderr,"loc 01 (%d,%d,%d)\n",tind,params.start_tind,params.Np);
      //fprintf(stderr,"loc 03\n");
      //	  fprintf(stderr,"EPV 01\n");
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        outputs.vertex_phasors.extractPhasorsVertices(
                ifx, inputs.E_s, inputs.H_s, tind, inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params);
      }
    }

    //fprintf(stderr,"Pos 02a:\n");
    if (inputs.params.source_mode == SourceMode::pulsed &&
        inputs.params.run_mode == RunMode::complete && inputs.params.exdetintegral) {
      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        spdlog::debug("Setting Ex_t, Ey_t");

        //First need to sum up the Ex and Ey values on a plane ready for FFT, remember that Ex_t and
        // Ey_t are in row-major format whilst Exy etc. are in column major format
        for (j = inputs.params.pml.Dyl; j < (J_tot - inputs.params.pml.Dyu); j++)
          for (i = inputs.params.pml.Dxl; i < (I_tot - inputs.params.pml.Dxu); i++) {
            int m = j - inputs.params.pml.Dyl +
                    (i - inputs.params.pml.Dxl) *
                            (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
            loop_variables.Ex_t.v[m][0] = inputs.E_s.xy[inputs.params.k_det_obs][j][i] +
                                          inputs.E_s.xz[inputs.params.k_det_obs][j][i];
            loop_variables.Ex_t.v[m][1] = 0.;
            loop_variables.Ey_t.v[m][0] = inputs.E_s.yx[inputs.params.k_det_obs][j][i] +
                                          inputs.E_s.yz[inputs.params.k_det_obs][j][i];
            loop_variables.Ey_t.v[m][1] = 0.;
          }

        fftw_execute(loop_variables.Ex_t.plan);
        fftw_execute(loop_variables.Ey_t.plan);

        //Iterate over each mode
        for (int im = 0; im < inputs.D_tilde.num_det_modes(); im++) {

          //Now go back to column-major
          for (j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl); j++)
            for (i = 0; i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
              int m = j + i * (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
              loop_variables.Ex_t.cm[j][i] =
                      loop_variables.Ex_t.v[m][0] + IMAGINARY_UNIT * loop_variables.Ex_t.v[m][1];
              loop_variables.Ey_t.cm[j][i] =
                      loop_variables.Ey_t.v[m][0] + IMAGINARY_UNIT * loop_variables.Ey_t.v[m][1];
            }

          //fprintf(stderr,"Pos 02a [3]:\n");
          //Now multiply the pupil, mostly the pupil is non-zero in only a elements
          for (j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl); j++)
            for (i = 0; i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
              loop_variables.Ex_t.cm[j][i] *= inputs.pupil[j][i] * inputs.D_tilde.x[j][i][im];
              loop_variables.Ey_t.cm[j][i] *= inputs.pupil[j][i] * inputs.D_tilde.y[j][i][im];
            }

            //now iterate over each frequency to extract phasors at
#pragma omp parallel default(shared) private(lambda_an_t, Idxt, Idyt, i, j, kprop, phaseTermE,     \
                                             cphaseTermE)
          {
#pragma omp for
            for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
              //wavelength in air
              lambda_an_t = LIGHT_V / inputs.f_ex_vec[ifx];
              //fprintf(stdout,"lambda_an_t = %e, LIGHT_V = %e, z_obs = %e\n",lambda_an_t,LIGHT_V,z_obs);
              Idxt = 0.;
              Idyt = 0.;

              //now loop over all angular frequencies
              for (j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl); j++)
                for (i = 0; i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
                  if ((lambda_an_t * inputs.f_vec.x[i] * lambda_an_t * inputs.f_vec.x[i] +
                       lambda_an_t * inputs.f_vec.y[j] * lambda_an_t * inputs.f_vec.y[j]) < 1) {

                    if (!inputs.params.air_interface_present) {
                      //This had to be fixed since we must take into account the refractive index of the medium.
                      kprop = exp(IMAGINARY_UNIT * inputs.params.z_obs * 2. * DCPI / lambda_an_t *
                                  loop_variables.refind *
                                  sqrt(1. -
                                       pow(lambda_an_t * inputs.f_vec.x[i] / loop_variables.refind,
                                           2.) -
                                       pow(lambda_an_t * inputs.f_vec.y[j] / loop_variables.refind,
                                           2.)));
                      //fprintf(stdout,"%d %d %e %e %e %e %e %e %e\n",i,j,f_vec.x[i],f_vec.y[j],real(kprop),imag(kprop),z_obs,DCPI,lambda_an_t);
                    } else {
                      kprop = exp(IMAGINARY_UNIT *
                                  (-inputs.params.air_interface + inputs.params.z_obs) * 2. * DCPI /
                                  lambda_an_t * loop_variables.refind *
                                  sqrt(1. -
                                       pow(lambda_an_t * inputs.f_vec.x[i] / loop_variables.refind,
                                           2.) -
                                       pow(lambda_an_t * inputs.f_vec.y[j] / loop_variables.refind,
                                           2.))) *
                              exp(IMAGINARY_UNIT * inputs.params.air_interface * 2. * DCPI /
                                  lambda_an_t *
                                  sqrt(1. - pow(lambda_an_t * inputs.f_vec.x[i], 2.) -
                                       pow(lambda_an_t * inputs.f_vec.y[j], 2.)));
                    }
                  } else
                    kprop = 0.;

                  Idxt += loop_variables.Ex_t.cm[j][i] * kprop;
                  Idyt += loop_variables.Ey_t.cm[j][i] * kprop;
                }
              phaseTermE =
                      fmod(inputs.f_ex_vec[ifx] * 2. * DCPI * ((double) tind) * inputs.params.dt,
                           2 * DCPI);
              cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. / ((double) inputs.params.Npe);

              outputs.ID.x[ifx][im] += Idxt * cphaseTermE;
              outputs.ID.y[ifx][im] += Idyt * cphaseTermE;

            }//end of loop on frequencies
          }  //end of pragma omp parallel
        }    //end of loop over each mode
      }
    }//end of section for calculating detector function

    //fprintf(stderr,"Pos 02b:\n");
    if (inputs.params.run_mode == RunMode::complete)
      if (inputs.params.dimension == THREE) {
        //extract the phasors just above the line
        FDTD.extract_phasors_in_plane(inputs.E_s, inputs.H_s, IJK_tot, inputs.K0.index + 1, tind,
                                      inputs.params);
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
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
    //fprintf(stderr,"Dimension = %d\n",params.dimension);
    /*
      for(k=0;k<(K_tot+1);k++)
      fprintf(stdout,"%e ",Exy[k][13][13]+Exz[k][13][13]);
      fprintf(stdout,"\n");
    */
    (void) n;// n is unused in FD derivatives – this silences the compiler warning

#pragma omp parallel default(shared) private(i, j, k, n, rho, k_loc, array_ind, Ca, Cb, Cc,        \
                                             alpha_l, beta_l, gamma_l, kappa_l, sigma_l, Enp1,     \
                                             Jnp1)//,ca_vec,cb_vec,eh_vec)
    {
      n = omp_get_thread_num();
      Enp1 = 0.0;
      array_ind = 0;

      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
          //FDTD, E_s.xy
#pragma omp for
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 1; j < J_tot; j++)
              for (i = 0; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (inputs.materials[k][j][i] || inputs.materials[k][j][i + 1]) {
                  //fprintf(stdout,"(%d,%d,%d,%d)\n",i,j,k,tind);
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
                             inputs.H_s.zy[k][j - 1][i] - inputs.H_s.zx[k][j - 1][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.xy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * loop_variables.J_s.xy[k][j][i] +
                                   beta_l * loop_variables.J_nm1.xy[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dy * loop_variables.J_c.xy[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.xy[k][j][i] +
                         beta_l * loop_variables.J_nm1.xy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.xy[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xy[k][j][i];

                  loop_variables.E_nm1.xy[k][j][i] = inputs.E_s.xy[k][j][i];
                  loop_variables.J_nm1.xy[k][j][i] = loop_variables.J_s.xy[k][j][i];
                  loop_variables.J_s.xy[k][j][i] = Jnp1;

                  //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,J_s.xy[k][j][i]);
                }

                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.xy[k][j][i] -= rho * (Enp1 + inputs.E_s.xy[k][j][i]);
                }

                inputs.E_s.xy[k][j][i] = Enp1;
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
                if (inputs.materials[k][j][i] || inputs.materials[k][j][i + 1]) {
                  //fprintf(stdout,"(%d,%d,%d,%d)\n",i,j,k,tind);
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
                //Enp1 = Ca*E_s.xy[k][j][i]+Cb*(H_s.zy[k][j][i] + H_s.zx[k][j][i] - H_s.zy[k][j-1][i] - H_s.zx[k][j-1][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.xy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * loop_variables.J_s.xy[k][j][i] +
                                   beta_l * loop_variables.J_nm1.xy[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dy * loop_variables.J_c.xy[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.xy[k][j][i] +
                         beta_l * loop_variables.J_nm1.xy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.xy[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xy[k][j][i];

                  loop_variables.E_nm1.xy[k][j][i] = inputs.E_s.xy[k][j][i];
                  loop_variables.J_nm1.xy[k][j][i] = loop_variables.J_s.xy[k][j][i];
                  loop_variables.J_s.xy[k][j][i] = Jnp1;

                  //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,J_s.xy[k][j][i]);
                }

                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.xy[k][j][i] -= rho * (Enp1 + inputs.E_s.xy[k][j][i]);
                }

                loop_variables.eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
                loop_variables.eh_vec[n][j][1] = 0.;
                PSTD.ca[n][j - 1] = Ca;
                PSTD.cb[n][j - 1] = Cb;
              }
              if (J_tot > 1) {
                j = 0;
                loop_variables.eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
                loop_variables.eh_vec[n][j][1] = 0.;
                first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_ey, PSTD.N_ey,
                                 inputs.E_s.xy.plan_f[n], inputs.E_s.xy.plan_b[n]);


                //fprintf(stdout,"(%d,%d) %d (of %d)\n",i,k,n,omp_get_num_threads());

                for (j = 1; j < J_tot; j++) {
                  inputs.E_s.xy[k][j][i] =
                          PSTD.ca[n][j - 1] * inputs.E_s.xy[k][j][i] +
                          PSTD.cb[n][j - 1] * loop_variables.eh_vec[n][j][0] / ((double) PSTD.N_ey);
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
            for (j = 0; j < loop_variables.J_tot_p1_bound; j++)
              for (i = 0; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                //use the average of material parameters between nodes
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
                /*if( materials[k][j][i] || materials[k][j][i+1])
      fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);
      if(tind==0)
      fprintf(stdout,"%d %d %e %e\n",i,k,Ca, Cb);*/
                Enp1 = Ca * inputs.E_s.xz[k][j][i] +
                       Cb * (inputs.H_s.yx[k - 1][j][i] + inputs.H_s.yz[k - 1][j][i] -
                             inputs.H_s.yx[k][j][i] - inputs.H_s.yz[k][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.xz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * loop_variables.J_s.xz[k][j][i] +
                                   beta_l * loop_variables.J_nm1.xz[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dz * loop_variables.J_c.xz[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.xz[k][j][i] +
                         beta_l * loop_variables.J_nm1.xz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.xz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xz[k][j][i];
                  loop_variables.E_nm1.xz[k][j][i] = inputs.E_s.xz[k][j][i];
                  loop_variables.J_nm1.xz[k][j][i] = loop_variables.J_s.xz[k][j][i];
                  loop_variables.J_s.xz[k][j][i] = Jnp1;
                }

                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.xz[k][j][i] -= rho * (Enp1 + inputs.E_s.xz[k][j][i]);
                }

                inputs.E_s.xz[k][j][i] = Enp1;
              }
          //FDTD, E_s.xz
        } else {
          //#pragma omp for
          for (j = 0; j < loop_variables.J_tot_p1_bound; j++)
#pragma omp for
            for (i = 0; i < I_tot; i++) {
              for (k = 1; k < K_tot; k++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                //use the average of material parameters between nodes
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
                /*if( materials[k][j][i] || materials[k][j][i+1])
      fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);
      if(tind==0)
      fprintf(stdout,"%d %d %e %e\n",i,k,Ca, Cb);*/
                //Enp1 = Ca*E_s.xz[k][j][i]+Cb*(H_s.yx[k-1][j][i] + H_s.yz[k-1][j][i] - H_s.yx[k][j][i] - H_s.yz[k][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.xz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * loop_variables.J_s.xz[k][j][i] +
                                   beta_l * loop_variables.J_nm1.xz[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dz * loop_variables.J_c.xz[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.xz[k][j][i] +
                         beta_l * loop_variables.J_nm1.xz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.xz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xz[k][j][i];
                  loop_variables.E_nm1.xz[k][j][i] = inputs.E_s.xz[k][j][i];
                  loop_variables.J_nm1.xz[k][j][i] = loop_variables.J_s.xz[k][j][i];
                  loop_variables.J_s.xz[k][j][i] = Jnp1;
                }

                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.xz[k][j][i] -= rho * (Enp1 + inputs.E_s.xz[k][j][i]);
                }

                loop_variables.eh_vec[n][k][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
                loop_variables.eh_vec[n][k][1] = 0.;
                PSTD.ca[n][k - 1] = Ca;
                PSTD.cb[n][k - 1] = Cb;
              }
              k = 0;
              loop_variables.eh_vec[n][k][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
              loop_variables.eh_vec[n][k][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_ez, PSTD.N_ez, inputs.E_s.xz.plan_f[n],
                               inputs.E_s.xz.plan_b[n]);

              for (k = 1; k < K_tot; k++) {
                inputs.E_s.xz[k][j][i] = PSTD.ca[n][k - 1] * inputs.E_s.xz[k][j][i] -
                                         PSTD.cb[n][k - 1] * loop_variables.eh_vec[n][k][0] / ((double) PSTD.N_ez);
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
            for (j = 0; j < loop_variables.J_tot_bound; j++)
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure) {
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.x[array_ind];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.x[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.x[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.x[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.x[array_ind];
                  kappa_l = inputs.matched_layer.kappa.x[array_ind];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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
                       Cb * (inputs.H_s.zx[k][j][i - 1] + inputs.H_s.zy[k][j][i - 1] -
                             inputs.H_s.zx[k][j][i] - inputs.H_s.zy[k][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * loop_variables.J_s.yx[k][j][i] +
                                   beta_l * loop_variables.J_nm1.yx[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx * loop_variables.J_c.yx[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yx[k][j][i] +
                         beta_l * loop_variables.J_nm1.yx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx[k][j][i];
                  loop_variables.E_nm1.yx[k][j][i] = inputs.E_s.yx[k][j][i];
                  loop_variables.J_nm1.yx[k][j][i] = loop_variables.J_s.yx[k][j][i];
                  loop_variables.J_s.yx[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yx[k][j][i] -= rho * (Enp1 + inputs.E_s.yx[k][j][i]);
                }

                inputs.E_s.yx[k][j][i] = Enp1;
              }
          //FDTD, E_s.yx
        } else {
#pragma omp for
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_tot_bound; j++) {
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure) {
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
                if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.x[array_ind];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.x[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.x[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.x[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.x[array_ind];
                  kappa_l = inputs.matched_layer.kappa.x[array_ind];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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


                //Enp1 = Ca*E_s.yx[k][j][i]+Cb*(H_s.zx[k][j][i-1] + H_s.zy[k][j][i-1] - H_s.zx[k][j][i] - H_s.zy[k][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * loop_variables.J_s.yx[k][j][i] +
                                   beta_l * loop_variables.J_nm1.yx[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx * loop_variables.J_c.yx[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yx[k][j][i] +
                         beta_l * loop_variables.J_nm1.yx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx[k][j][i];
                  loop_variables.E_nm1.yx[k][j][i] = inputs.E_s.yx[k][j][i];
                  loop_variables.J_nm1.yx[k][j][i] = loop_variables.J_s.yx[k][j][i];
                  loop_variables.J_s.yx[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yx[k][j][i] -= rho * (Enp1 + inputs.E_s.yx[k][j][i]);
                }

                loop_variables.eh_vec[n][i][0] = inputs.H_s.zx[k][j][i] + inputs.H_s.zy[k][j][i];
                loop_variables.eh_vec[n][i][1] = 0.;
                PSTD.ca[n][i - 1] = Ca;
                PSTD.cb[n][i - 1] = Cb;
              }
              i = 0;
              loop_variables.eh_vec[n][i][0] = inputs.H_s.zx[k][j][i] + inputs.H_s.zy[k][j][i];
              loop_variables.eh_vec[n][i][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_ex, PSTD.N_ex, inputs.E_s.yx.plan_f[n],
                               inputs.E_s.yx.plan_b[n]);

              for (i = 1; i < I_tot; i++) {
                inputs.E_s.yx[k][j][i] = PSTD.ca[n][i - 1] * inputs.E_s.yx[k][j][i] -
                                         PSTD.cb[n][i - 1] * loop_variables.eh_vec[n][i][0] / ((double) PSTD.N_ex);
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
            for (j = 0; j < loop_variables.J_tot_bound; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      Ca = Ca + inputs.Cmaterial.a.z[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.z[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.z[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.z[k_loc];
                  kappa_l = inputs.matched_layer.kappa.z[k_loc];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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

                //fprintf(stderr,"[%d %d %d]Ca: %e, Cb: %e, Cc: %e, alpha: %e, beta: %e, gamme: %e\n",i,j,k,Ca,Cb,Cc,alpha_l,beta_l,gamma_l);
                Enp1 = Ca * inputs.E_s.yz[k][j][i] +
                       Cb * (inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i] -
                             inputs.H_s.xy[k - 1][j][i] - inputs.H_s.xz[k - 1][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * loop_variables.J_s.yz[k][j][i] +
                                   beta_l * loop_variables.J_nm1.yz[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dz * loop_variables.J_c.yz[k][j][i];

                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yz[k][j][i] +
                         beta_l * loop_variables.J_nm1.yz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz[k][j][i];
                  loop_variables.E_nm1.yz[k][j][i] = inputs.E_s.yz[k][j][i];
                  loop_variables.J_nm1.yz[k][j][i] = loop_variables.J_s.yz[k][j][i];
                  loop_variables.J_s.yz[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yz[k][j][i] -= rho * (Enp1 + inputs.E_s.yz[k][j][i]);
                }

                inputs.E_s.yz[k][j][i] = Enp1;
              }
          //FDTD, E_s.yz
        } else {
#pragma omp for
          for (j = 0; j < loop_variables.J_tot_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (k = 1; k < K_tot; k++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      Ca = Ca + inputs.Cmaterial.a.z[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.z[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.z[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
                  sigma_l = inputs.matched_layer.sigma.z[k_loc];
                  kappa_l = inputs.matched_layer.kappa.z[k_loc];
                  alpha_l = inputs.matched_layer.alpha[k_loc];
                  beta_l = inputs.matched_layer.beta[k_loc];
                  gamma_l = inputs.matched_layer.gamma[k_loc];
                  if (inputs.materials[k][j][i] || inputs.materials[k][min(J_tot, j + 1)][i]) {
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
                      alpha_l += inputs.alpha[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      beta_l += inputs.beta[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
                      gamma_l += inputs.gamma[inputs.materials[k][min(J_tot, j + 1)][i] - 1];
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

                //fprintf(stderr,"[%d %d %d]Ca: %e, Cb: %e, Cc: %e, alpha: %e, beta: %e, gamme: %e\n",i,j,k,Ca,Cb,Cc,alpha_l,beta_l,gamma_l);
                //Enp1 = Ca*E_s.yz[k][j][i]+Cb*(H_s.xy[k][j][i] + H_s.xz[k][j][i] - H_s.xy[k-1][j][i] - H_s.xz[k-1][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.yz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * loop_variables.J_s.yz[k][j][i] +
                                   beta_l * loop_variables.J_nm1.yz[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dz * loop_variables.J_c.yz[k][j][i];

                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.yz[k][j][i] +
                         beta_l * loop_variables.J_nm1.yz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.yz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz[k][j][i];
                  loop_variables.E_nm1.yz[k][j][i] = inputs.E_s.yz[k][j][i];
                  loop_variables.J_nm1.yz[k][j][i] = loop_variables.J_s.yz[k][j][i];
                  loop_variables.J_s.yz[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.yz[k][j][i] -= rho * (Enp1 + inputs.E_s.yz[k][j][i]);
                }

                loop_variables.eh_vec[n][k][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
                loop_variables.eh_vec[n][k][1] = 0.;
                PSTD.ca[n][k - 1] = Ca;
                PSTD.cb[n][k - 1] = Cb;
              }
              k = 0;
              loop_variables.eh_vec[n][k][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
              loop_variables.eh_vec[n][k][1] = 0.;
              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_ez, PSTD.N_ez, inputs.E_s.yz.plan_f[n],
                               inputs.E_s.yz.plan_b[n]);


              for (k = 1; k < K_tot; k++) {
                inputs.E_s.yz[k][j][i] = PSTD.ca[n][k - 1] * inputs.E_s.yz[k][j][i] +
                                         PSTD.cb[n][k - 1] * loop_variables.eh_vec[n][k][0] / ((double) PSTD.N_ez);
                //E_s.yz[k][j][i] = Enp1;
              }
            }
          //PSTD, E_s.yz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
      }  //if(params.dimension==THREE || params.dimension==TE)

      //fprintf(stderr,"Pos 07:\n");
      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
#pragma omp for
          //E_s.zx updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_tot_p1_bound; j++)
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
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
                      if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.x[array_ind];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.x[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.x[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.x[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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

                /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);*/
                Enp1 = Ca * inputs.E_s.zx[k][j][i] +
                       Cb * (inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i] -
                             inputs.H_s.yx[k][j][i - 1] - inputs.H_s.yz[k][j][i - 1]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * loop_variables.J_s.zx[k][j][i] +
                                   beta_l * loop_variables.J_nm1.zx[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx * loop_variables.J_c.zx[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zx[k][j][i] +
                         beta_l * loop_variables.J_nm1.zx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
                  loop_variables.E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
                  loop_variables.J_nm1.zx[k][j][i] = loop_variables.J_s.zx[k][j][i];
                  loop_variables.J_s.zx[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]);
                }

                inputs.E_s.zx[k][j][i] = Enp1;
              }
          //FDTD, E_s.zx
        } else {
#pragma omp for
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_tot_p1_bound; j++) {
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;

                //use the average of material parameters between nodes
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
                      if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.x[array_ind];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.x[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.x[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.x[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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

                /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_conductive:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_conductive,rho,is_disp,params.is_disp_ml);*/
                //Enp1 = Ca*E_s.zx[k][j][i]+Cb*(H_s.yx[k][j][i] + H_s.yz[k][j][i] - H_s.yx[k][j][i-1] - H_s.yz[k][j][i-1]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * loop_variables.J_s.zx[k][j][i] +
                                   beta_l * loop_variables.J_nm1.zx[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dx * loop_variables.J_c.zx[k][j][i];
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zx[k][j][i] +
                         beta_l * loop_variables.J_nm1.zx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
                  loop_variables.E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
                  loop_variables.J_nm1.zx[k][j][i] = loop_variables.J_s.zx[k][j][i];
                  loop_variables.J_s.zx[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]);
                }

                loop_variables.eh_vec[n][i][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
                loop_variables.eh_vec[n][i][1] = 0.;
                PSTD.ca[n][i - 1] = Ca;
                PSTD.cb[n][i - 1] = Cb;
              }
              i = 0;
              loop_variables.eh_vec[n][i][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
              loop_variables.eh_vec[n][i][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_ex, PSTD.N_ex, inputs.E_s.zx.plan_f[n],
                               inputs.E_s.zx.plan_b[n]);

              for (i = 1; i < I_tot; i++) {
                inputs.E_s.zx[k][j][i] = PSTD.ca[n][i - 1] * inputs.E_s.zx[k][j][i] +
                                         PSTD.cb[n][i - 1] * loop_variables.eh_vec[n][i][0] / ((double) PSTD.N_ex);
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
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                  if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = i;
              else
                array_ind = (I_tot + 1) * k_loc + i;

              //use the average of material parameters between nodes
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


              if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
              if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                Enp1 += Cc * loop_variables.E_nm1.zx[k][j][i] -
                        1. / 2. * Cb * inputs.params.delta.dx *
                                ((1 + alpha_l) * loop_variables.J_s.zx[k][j][i] +
                                 beta_l * loop_variables.J_nm1.zx[k][j][i]);
              if (loop_variables.is_conductive && rho)
                Enp1 += Cb * inputs.params.delta.dx * loop_variables.J_c.zx[k][j][i];

              if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * loop_variables.J_s.zx[k][j][i] +
                       beta_l * loop_variables.J_nm1.zx[k][j][i] +
                       kappa_l * gamma_l / (2. * inputs.params.dt) *
                               (Enp1 - loop_variables.E_nm1.zx[k][j][i]);
                Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
                loop_variables.E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
                loop_variables.J_nm1.zx[k][j][i] = loop_variables.J_s.zx[k][j][i];
                loop_variables.J_s.zx[k][j][i] = Jnp1;
              }
              if (loop_variables.is_conductive && rho) {
                loop_variables.J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]);
              }

              inputs.E_s.zx[k][j][i] = Enp1;
            }
      }
      //fprintf(stderr,"Pos 08:\n");
      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
          //FDTD, E_s.zy
#pragma omp for
          //E_s.zy updates
          for (k = 0; k < K_tot; k++)
            for (j = 1; j < J_tot; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
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
                      if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.y[array_ind];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.y[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.y[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.y[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
                       Cb * (inputs.H_s.xy[k][j - 1][i] + inputs.H_s.xz[k][j - 1][i] -
                             inputs.H_s.xy[k][j][i] - inputs.H_s.xz[k][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * loop_variables.J_s.zy[k][j][i] +
                                   beta_l * loop_variables.J_nm1.zy[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dy * loop_variables.J_c.zy[k][j][i];

                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zy[k][j][i] +
                         beta_l * loop_variables.J_nm1.zy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zy[k][j][i]);

                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
                  loop_variables.E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
                  loop_variables.J_nm1.zy[k][j][i] = loop_variables.J_s.zy[k][j][i];
                  loop_variables.J_s.zy[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]);
                }
                inputs.E_s.zy[k][j][i] = Enp1;
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;

                //use the average of material parameters between nodes
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
                      if (inputs.params.is_disp_ml) Cc = Cc + inputs.C.c.y[array_ind];
                    } else {
                      Ca = Ca + inputs.Cmaterial.a.y[inputs.materials[k + 1][j][i] - 1];
                      Cb = Cb + inputs.Cmaterial.b.y[inputs.materials[k + 1][j][i] - 1];
                      Cc = Cc + inputs.Cmaterial.c.y[inputs.materials[k + 1][j][i] - 1];
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
                  if (loop_variables.is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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


                //Enp1 = Ca*E_s.zy[k][j][i]+Cb*(H_s.xy[k][j-1][i] + H_s.xz[k][j-1][i] - H_s.xy[k][j][i] - H_s.xz[k][j][i]);
                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * loop_variables.E_nm1.zy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * loop_variables.J_s.zy[k][j][i] +
                                   beta_l * loop_variables.J_nm1.zy[k][j][i]);
                if (loop_variables.is_conductive && rho)
                  Enp1 += Cb * inputs.params.delta.dy * loop_variables.J_c.zy[k][j][i];

                if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * loop_variables.J_s.zy[k][j][i] +
                         beta_l * loop_variables.J_nm1.zy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) *
                                 (Enp1 - loop_variables.E_nm1.zy[k][j][i]);

                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
                  loop_variables.E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
                  loop_variables.J_nm1.zy[k][j][i] = loop_variables.J_s.zy[k][j][i];
                  loop_variables.J_s.zy[k][j][i] = Jnp1;
                }
                if (loop_variables.is_conductive && rho) {
                  loop_variables.J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]);
                }

                loop_variables.eh_vec[n][j][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
                loop_variables.eh_vec[n][j][1] = 0.;
                PSTD.ca[n][j - 1] = Ca;
                PSTD.cb[n][j - 1] = Cb;
              }
              if (J_tot > 1) {
                j = 0;
                loop_variables.eh_vec[n][j][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
                loop_variables.eh_vec[n][j][1] = 0.;
                first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_ey, PSTD.N_ey,
                                 inputs.E_s.zy.plan_f[n], inputs.E_s.zy.plan_b[n]);
              }
              for (j = 1; j < J_tot; j++) {
                inputs.E_s.zy[k][j][i] = PSTD.ca[n][j - 1] * inputs.E_s.zy[k][j][i] -
                                         PSTD.cb[n][j - 1] * loop_variables.eh_vec[n][j][0] / ((double) PSTD.N_ey);
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
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                  if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = j;
              else
                array_ind = (J_tot + 1) * k_loc + j;

              //use the average of material parameters between nodes
              if (!inputs.materials[k][j][i]) {
                Ca = inputs.C.a.y[array_ind];
                Cb = inputs.C.b.y[array_ind];
                if (inputs.params.is_disp_ml) Cc = inputs.C.c.y[array_ind];
                else
                  Cc = 0.;
                if (loop_variables.is_conductive) rho = inputs.rho_cond.y[array_ind];
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

              if (loop_variables.is_disp || inputs.params.is_disp_ml) {
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
              if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l)
                Enp1 += Cc * loop_variables.E_nm1.zy[k][j][i] -
                        1. / 2. * Cb * inputs.params.delta.dy *
                                ((1 + alpha_l) * loop_variables.J_s.zy[k][j][i] +
                                 beta_l * loop_variables.J_nm1.zy[k][j][i]);
              if (loop_variables.is_conductive && rho)
                Enp1 += Cb * inputs.params.delta.dy * loop_variables.J_c.zy[k][j][i];

              if ((loop_variables.is_disp || inputs.params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * loop_variables.J_s.zy[k][j][i] +
                       beta_l * loop_variables.J_nm1.zy[k][j][i] +
                       kappa_l * gamma_l / (2. * inputs.params.dt) *
                               (Enp1 - loop_variables.E_nm1.zy[k][j][i]);

                Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
                loop_variables.E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
                loop_variables.J_nm1.zy[k][j][i] = loop_variables.J_s.zy[k][j][i];
                loop_variables.J_s.zy[k][j][i] = Jnp1;
              }
              if (loop_variables.is_conductive && rho) {
                loop_variables.J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]);
              }

              inputs.E_s.zy[k][j][i] = Enp1;
            }
      }
    }//end of parallel section
    //fprintf(stderr,"Pos 09:\n");
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
    /********************/

    //update terms for self consistency across scattered/total interface - E updates##
    if (inputs.params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase =
              exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_H, 2. * DCPI));
      double commonAmplitude = linear_ramp(time_H);
      for (k = (inputs.K0.index); k <= (inputs.K1.index); k++)
        for (j = (inputs.J0.index); j <= (inputs.J1.index); j++) {
          if (inputs.I0.apply) {//Perform across I0

            if (!inputs.params.is_multilayer) array_ind = inputs.I0.index;
            else
              array_ind = (I_tot + 1) * k + inputs.I0.index;

            if (k < (inputs.K1.index) ||
                inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
              inputs.E_s.zx[k][j][inputs.I0.index] =
                      inputs.E_s.zx[k][j][inputs.I0.index] -
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][2] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][2]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.zx[k][j][inputs.I0.index] +=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][2] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][2]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.zx[k][j][inputs.I0.index] +=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][2] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][2]));
            }
            if (j < (inputs.J1.index)) {
              inputs.E_s.yx[k][j][inputs.I0.index] =
                      inputs.E_s.yx[k][j][inputs.I0.index] +
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][3] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][3]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.yx[k][j][inputs.I0.index] -=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][3] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][3]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.yx[k][j][inputs.I0.index] -=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][3] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][3]));
            }
          }
          if (inputs.I1.apply) {//Perform across I1

            if (!inputs.params.is_multilayer) array_ind = inputs.I1.index;
            else
              array_ind = (I_tot + 1) * k + inputs.I1.index;

            if (k < (inputs.K1.index) ||
                inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
              inputs.E_s.zx[k][j][inputs.I1.index] =
                      inputs.E_s.zx[k][j][inputs.I1.index] +
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][6] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][6]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.zx[k][j][inputs.I1.index] -=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][6] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][6]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.zx[k][j][inputs.I1.index] -=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][6] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][6]));
            }
            if (j < (inputs.J1.index)) {
              inputs.E_s.yx[k][j][inputs.I1.index] =
                      inputs.E_s.yx[k][j][inputs.I1.index] -
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][7] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][7]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.yx[k][j][inputs.I1.index] +=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][7] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][7]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.yx[k][j][inputs.I1.index] +=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][7] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)]
                                                                  [j - (inputs.J0.index)][7]));
            }
          }
        }

      for (k = (inputs.K0.index); k <= (inputs.K1.index); k++)
        for (i = (inputs.I0.index); i <= (inputs.I1.index); i++) {
          if (inputs.J0.apply) {//Perform across J0
            if (k < (inputs.K1.index) ||
                inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {

              if (!inputs.params.is_multilayer) array_ind = inputs.J0.index;
              else
                array_ind = (J_tot + 1) * k + inputs.J0.index;

              inputs.E_s.zy[k][(inputs.J0.index)][i] =
                      inputs.E_s.zy[k][(inputs.J0.index)][i] +
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][2] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][2]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.zy[k][(inputs.J0.index)][i] -=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][2]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.zy[k][(inputs.J0.index)][i] -=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][2]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xy[k][(inputs.J0.index)][i] =
                      inputs.E_s.xy[k][(inputs.J0.index)][i] -
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][3] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][3]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.xy[k][(inputs.J0.index)][i] +=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][3]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.xy[k][(inputs.J0.index)][i] +=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][3]));
            }
          }
          if (inputs.J1.apply) {//Perform across J1

            if (!inputs.params.is_multilayer) array_ind = inputs.J1.index;
            else
              array_ind = (J_tot + 1) * k + inputs.J1.index;

            if (k < (inputs.K1.index) ||
                inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
              inputs.E_s.zy[k][(inputs.J1.index)][i] =
                      inputs.E_s.zy[k][(inputs.J1.index)][i] -
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][6] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][6]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.zy[k][(inputs.J1.index)][i] +=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][6]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.zy[k][(inputs.J1.index)][i] -=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][6]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xy[k][(inputs.J1.index)][i] =
                      inputs.E_s.xy[k][(inputs.J1.index)][i] +
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][7] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][7]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.xy[k][(inputs.J1.index)][i] -=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][7]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.xy[k][(inputs.J1.index)][i] +=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] /
                        (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)]
                                                                  [i - (inputs.I0.index)][7]));
            }
          }
        }

      for (j = (inputs.J0.index); j <= (inputs.J1.index); j++)
        for (i = (inputs.I0.index); i <= (inputs.I1.index); i++) {
          if (inputs.K0.apply) {//Perform across K0
            if (j < (inputs.J1.index)) {
              inputs.E_s.yz[(inputs.K0.index)][j][i] =
                      inputs.E_s.yz[(inputs.K0.index)][j][i] -
                      inputs.C.b.z[inputs.K0.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][2]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.yz[(inputs.K0.index)][j][i] +=
                        inputs.rho_cond.z[(inputs.K0.index)] * inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][2]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.yz[(inputs.K0.index)][j][i] -=
                        inputs.matched_layer.kappa.z[(inputs.K0.index)] *
                        inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                        inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][2]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xz[(inputs.K0.index)][j][i] =
                      inputs.E_s.xz[(inputs.K0.index)][j][i] +
                      inputs.C.b.z[inputs.K0.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][3]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.xz[(inputs.K0.index)][j][i] -=
                        inputs.rho_cond.z[(inputs.K0.index)] * inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][3]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.xz[(inputs.K0.index)][j][i] +=
                        inputs.matched_layer.kappa.z[(inputs.K0.index)] *
                        inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                        inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][3]));
            }
          }
          if (inputs.K1.apply) {//Perform across K1
            if (j < (inputs.J1.index)) {
              inputs.E_s.yz[(inputs.K1.index)][j][i] =
                      inputs.E_s.yz[(inputs.K1.index)][j][i] +
                      inputs.C.b.z[inputs.K1.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][6] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][6]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.yz[(inputs.K1.index)][j][i] -=
                        inputs.rho_cond.z[(inputs.K1.index)] * inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][6]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.yz[(inputs.K1.index)][j][i] +=
                        inputs.matched_layer.kappa.z[(inputs.K1.index)] *
                        inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                        inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][6]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xz[(inputs.K1.index)][j][i] =
                      inputs.E_s.xz[(inputs.K1.index)][j][i] -
                      inputs.C.b.z[inputs.K1.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][7] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][7]));
              if (loop_variables.is_conductive)
                loop_variables.J_c.xz[(inputs.K1.index)][j][i] +=
                        inputs.rho_cond.z[(inputs.K1.index)] * inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][7]));
              if (inputs.params.is_disp_ml)
                loop_variables.J_s.xz[(inputs.K1.index)][j][i] -=
                        inputs.matched_layer.kappa.z[(inputs.K1.index)] *
                        inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) *
                        inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                  [i - (inputs.I0.index)][7]));
            }
          }
        }
      outputs.H.ft = real(commonAmplitude * commonPhase);
    } else if (inputs.params.source_mode == SourceMode::pulsed) {//pulsed

      if (J_tot == 0) {
        j = 0;
        for (i = 0; i < (I_tot + 1); i++) {
          inputs.E_s.yz[inputs.K0.index][j][i] =
                  inputs.E_s.yz[inputs.K0.index][j][i] -
                  inputs.C.b.z[inputs.K0.index] *
                          real((inputs.Ksource.real[0][i - (inputs.I0.index)][2] +
                                IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][2]) *
                               (-1.0 * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT *
                                   fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                        2. * DCPI))) *
                          exp(-1.0 * DCPI *
                              pow((time_H - inputs.params.to_l +
                                   inputs.params.delta.dz / LIGHT_V / 2.) /
                                          (inputs.params.hwhm),
                                  2));
          //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
          if (loop_variables.is_conductive)
            loop_variables.J_c.yz[inputs.K0.index][j][i] +=
                    inputs.rho_cond.z[inputs.K0.index] * inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[0][i - (inputs.I0.index)][2] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][2]) *
                         (-1.0 * IMAGINARY_UNIT) *
                         exp(-IMAGINARY_UNIT *
                             fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                  2. * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) /
                                    (inputs.params.hwhm),
                            2));
          //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
          if (inputs.params.is_disp_ml) {
            loop_variables.J_s.yz[inputs.K0.index][j][i] -=
                    inputs.matched_layer.kappa.z[inputs.K0.index] *
                    inputs.matched_layer.gamma[inputs.K0.index] / (2. * inputs.params.dt) *
                    inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[0][i - (inputs.I0.index)][2] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][2]) *
                         (-1.0 * IMAGINARY_UNIT) *
                         exp(-IMAGINARY_UNIT *
                             fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                  2. * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) /
                                    (inputs.params.hwhm),
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
            inputs.E_s.yz[inputs.K0.index][j][i] =
                    inputs.E_s.yz[inputs.K0.index][j][i] -
                    inputs.C.b.z[inputs.K0.index] *
                            real((inputs.Ksource
                                          .real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                                  IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                      [i - (inputs.I0.index)][2]) *
                                 (-1.0 * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT *
                                     fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                          2. * DCPI))) *
                            exp(-1.0 * DCPI *
                                pow((time_H - inputs.params.to_l +
                                     inputs.params.delta.dz / LIGHT_V / 2.) /
                                            (inputs.params.hwhm),
                                    2));
            //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (loop_variables.is_conductive)
              loop_variables.J_c.yz[inputs.K0.index][j][i] +=
                      inputs.rho_cond.z[inputs.K0.index] * inputs.C.b.z[inputs.K0.index] *
                      real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                            IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                [i - (inputs.I0.index)][2]) *
                           (-1.0 * IMAGINARY_UNIT) *
                           exp(-IMAGINARY_UNIT *
                               fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                    2. * DCPI))) *
                      exp(-1.0 * DCPI *
                          pow((time_H - inputs.params.to_l +
                               inputs.params.delta.dz / LIGHT_V / 2.) /
                                      (inputs.params.hwhm),
                              2));
            //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (inputs.params.is_disp_ml) {
              loop_variables.J_s.yz[inputs.K0.index][j][i] -=
                      inputs.matched_layer.kappa.z[inputs.K0.index] *
                      inputs.matched_layer.gamma[inputs.K0.index] / (2. * inputs.params.dt) *
                      inputs.C.b.z[inputs.K0.index] *
                      real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                            IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                [i - (inputs.I0.index)][2]) *
                           (-1.0 * IMAGINARY_UNIT) *
                           exp(-IMAGINARY_UNIT *
                               fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                    2. * DCPI))) *
                      exp(-1.0 * DCPI *
                          pow((time_H - inputs.params.to_l +
                               inputs.params.delta.dz / LIGHT_V / 2.) /
                                      (inputs.params.hwhm),
                              2));
              //J_s.yz[(int)K0[0]][j][i] -= matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            }
          }
      for (j = 0; j < (J_tot + 1); j++)
        for (i = 0; i < I_tot; i++) {
          inputs.E_s.xz[inputs.K0.index][j][i] =
                  inputs.E_s.xz[inputs.K0.index][j][i] +
                  inputs.C.b.z[inputs.K0.index] *
                          real((inputs.Ksource
                                        .real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                                IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                    [i - (inputs.I0.index)][3]) *
                               (-1.0 * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT *
                                   fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                        2 * DCPI))) *
                          exp(-1.0 * DCPI *
                              pow((time_H - inputs.params.to_l +
                                   inputs.params.delta.dz / LIGHT_V / 2.) /
                                          (inputs.params.hwhm),
                                  2));
          //E_s.xz[(int)K0[0]][j][i] = E_s.xz[(int)K0[0]][j][i] + C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
          if (loop_variables.is_conductive)
            loop_variables.J_c.xz[inputs.K0.index][j][i] -=
                    inputs.rho_cond.z[inputs.K0.index] * inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                              [i - (inputs.I0.index)][3]) *
                         (-1.0 * IMAGINARY_UNIT) *
                         exp(-IMAGINARY_UNIT *
                             fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                  2 * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) /
                                    (inputs.params.hwhm),
                            2));
          //J_c.xz[(int)K0[0]][j][i] -= rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
          if (inputs.params.is_disp_ml)
            loop_variables.J_s.xz[inputs.K0.index][j][i] +=
                    inputs.matched_layer.kappa.z[inputs.K0.index] *
                    inputs.matched_layer.gamma[inputs.K0.index] / (2. * inputs.params.dt) *
                    inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                              [i - (inputs.I0.index)][3]) *
                         (-1.0 * IMAGINARY_UNIT) *
                         exp(-IMAGINARY_UNIT *
                             fmod(inputs.params.omega_an * (time_H - inputs.params.to_l),
                                  2 * DCPI))) *
                    exp(-1.0 * DCPI *
                        pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) /
                                    (inputs.params.hwhm),
                            2));
          //J_s.xz[(int)K0[0]][j][i] += matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
        }
      //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
      outputs.H.ft =
              real((-1.0 * IMAGINARY_UNIT) *
                   exp(-IMAGINARY_UNIT *
                       fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
              exp(-1.0 * DCPI *
                  pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) /
                              (inputs.params.hwhm),
                      2));
      //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
    }
    //fprintf(stderr,"Pos 10:\n");

    //end of source terms
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    /********************/
    //begin parallel
#pragma omp parallel default(shared) private(i, j, k, n, k_loc, array_ind)//,ca_vec,cb_vec,eh_vec)
    {
      n = omp_get_thread_num();

      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
//FDTD, H_s.xz
#pragma omp for
          //H_s.xz updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_tot_bound; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }

                if (!inputs.materials[k][j][i])
                  inputs.H_s.xz[k][j][i] =
                          inputs.D.a.z[k_loc] * inputs.H_s.xz[k][j][i] +
                          inputs.D.b.z[k_loc] *
                                  (inputs.E_s.yx[k + 1][j][i] + inputs.E_s.yz[k + 1][j][i] -
                                   inputs.E_s.yx[k][j][i] - inputs.E_s.yz[k][j][i]);
                else
                  inputs.H_s.xz[k][j][i] =
                          inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.xz[k][j][i] +
                          inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.yx[k + 1][j][i] + inputs.E_s.yz[k + 1][j][i] -
                                   inputs.E_s.yx[k][j][i] - inputs.E_s.yz[k][j][i]);
              }
          //FDTD, H_s.xz
        } else {
#pragma omp for
          //H_s.xz updates
          for (j = 0; j < loop_variables.J_tot_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }

                if (!inputs.materials[k][j][i]) {
                  PSTD.ca[n][k] = inputs.D.a.z[k_loc];
                  PSTD.cb[n][k] = inputs.D.b.z[k_loc];
                  //H_s.xz[k][j][i] = D.a.z[k_loc]*H_s.xz[k][j][i]+D.b.z[k_loc]*(E_s.yx[k+1][j][i] + E_s.yz[k+1][j][i] - E_s.yx[k][j][i] - E_s.yz[k][j][i]);
                } else {
                  PSTD.ca[n][k] = inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1];
                  PSTD.cb[n][k] = inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1];
                  //H_s.xz[k][j][i] = Dmaterial.Da.z[materials[k][j][i]-1]*H_s.xz[k][j][i]+Dmaterial.Db.z[materials[k][j][i]-1]*(E_s.yx[k+1][j][i] + E_s.yz[k+1][j][i] - E_s.yx[k][j][i] - E_s.yz[k][j][i]);
                }

                loop_variables.eh_vec[n][k][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
                loop_variables.eh_vec[n][k][1] = 0.;
              }
              k = K_tot;
              loop_variables.eh_vec[n][k][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
              loop_variables.eh_vec[n][k][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_hz, PSTD.N_hz, inputs.H_s.xz.plan_f[n],
                               inputs.H_s.xz.plan_b[n]);

              for (k = 0; k < K_tot; k++) {
                inputs.H_s.xz[k][j][i] = PSTD.ca[n][k] * inputs.H_s.xz[k][j][i] +
                                         PSTD.cb[n][k] * loop_variables.eh_vec[n][k][0] / ((double) PSTD.N_hz);
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.xy[k][j][i] =
                          inputs.D.a.y[array_ind] * inputs.H_s.xy[k][j][i] +
                          inputs.D.b.y[array_ind] *
                                  (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] -
                                   inputs.E_s.zy[k][j + 1][i] - inputs.E_s.zx[k][j + 1][i]);
                else
                  inputs.H_s.xy[k][j][i] =
                          inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.xy[k][j][i] +
                          inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] -
                                   inputs.E_s.zy[k][j + 1][i] - inputs.E_s.zx[k][j + 1][i]);
              }
          //FDTD, H_s.xy
        } else {
#pragma omp for
          //H_s.xy updates
          for (k = 0; k < K_tot; k++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (j = 0; j < J_tot; j++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca[n][j] = inputs.D.a.y[array_ind];
                  PSTD.cb[n][j] = inputs.D.b.y[array_ind];
                  //		H_s.xy[k][j][i] = D.a.y[array_ind]*H_s.xy[k][j][i]+D.b.y[array_ind]*(E_s.zy[k][j][i] + E_s.zx[k][j][i] - E_s.zy[k][j+1][i] - E_s.zx[k][j+1][i]);
                } else {
                  PSTD.ca[n][j] = inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1];
                  PSTD.cb[n][j] = inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1];
                  //		H_s.xy[k][j][i] = Dmaterial.Da.y[materials[k][j][i]-1]*H_s.xy[k][j][i]+Dmaterial.Db.y[materials[k][j][i]-1]*(E_s.zy[k][j][i] + E_s.zx[k][j][i] - E_s.zy[k][j+1][i] - E_s.zx[k][j+1][i]);
                }

                loop_variables.eh_vec[n][j][0] = inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i];
                loop_variables.eh_vec[n][j][1] = 0.;
              }
              j = J_tot;
              loop_variables.eh_vec[n][j][0] = inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i];
              loop_variables.eh_vec[n][j][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_hy, PSTD.N_hy, inputs.H_s.xy.plan_f[n],
                               inputs.H_s.xy.plan_b[n]);

              for (j = 0; j < J_tot; j++) {
                inputs.H_s.xy[k][j][i] = PSTD.ca[n][j] * inputs.H_s.xy[k][j][i] -
                                         PSTD.cb[n][j] * loop_variables.eh_vec[n][j][0] / ((double) PSTD.N_hy);
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
    fprintf(stdout,"%e ",eh_vec[n][j][0]/((double) PSTD.N_ey));
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
            for (j = 0; j < loop_variables.J_tot_p1_bound; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.yx[k][j][i] =
                          inputs.D.a.x[array_ind] * inputs.H_s.yx[k][j][i] +
                          inputs.D.b.x[array_ind] *
                                  (inputs.E_s.zx[k][j][i + 1] + inputs.E_s.zy[k][j][i + 1] -
                                   inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j][i]);
                else {
                  inputs.H_s.yx[k][j][i] =
                          inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.yx[k][j][i] +
                          inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.zx[k][j][i + 1] + inputs.E_s.zy[k][j][i + 1] -
                                   inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j][i]);
                }
              }
          //FDTD, H_s.yx
        } else {
#pragma omp for
          //H_s.yx updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < loop_variables.J_tot_p1_bound; j++) {
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca[n][i] = inputs.D.a.x[array_ind];
                  PSTD.cb[n][i] = inputs.D.b.x[array_ind];
                  //		H_s.yx[k][j][i] = D.a.x[array_ind]*H_s.yx[k][j][i]+D.b.x[array_ind]*(E_s.zx[k][j][i+1] + E_s.zy[k][j][i+1] - E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                } else {
                  PSTD.ca[n][i] = inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1];
                  PSTD.cb[n][i] = inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1];
                  //	H_s.yx[k][j][i] = Dmaterial.Da.x[materials[k][j][i]-1]*H_s.yx[k][j][i]+Dmaterial.Db.x[materials[k][j][i]-1]*(E_s.zx[k][j][i+1] + E_s.zy[k][j][i+1] - E_s.zx[k][j][i] - E_s.zy[k][j][i]);
                }

                loop_variables.eh_vec[n][i][0] = inputs.E_s.zx[k][j][i] + inputs.E_s.zy[k][j][i];
                loop_variables.eh_vec[n][i][1] = 0.;
              }
              i = I_tot;
              loop_variables.eh_vec[n][i][0] = inputs.E_s.zx[k][j][i] + inputs.E_s.zy[k][j][i];
              loop_variables.eh_vec[n][i][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_hx, PSTD.N_hx, inputs.H_s.yx.plan_f[n],
                               inputs.H_s.yx.plan_b[n]);

              for (i = 0; i < I_tot; i++) {
                inputs.H_s.yx[k][j][i] = PSTD.ca[n][i] * inputs.H_s.yx[k][j][i] +
                                         PSTD.cb[n][i] * loop_variables.eh_vec[n][i][0] / ((double) PSTD.N_hx);
              }
            }
          //PSTD, H_s.yx
        }

        if (solver_method == SolverMethod::FiniteDifference) {
//FDTD, H_s.yz
#pragma omp for
          //H_s.yz updates
          for (k = 0; k < K_tot; k++) {
            for (j = 0; j < loop_variables.J_tot_p1_bound; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.materials[k][j][i]) {
                  /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                  inputs.H_s.yz[k][j][i] =
                          inputs.D.a.z[k_loc] * inputs.H_s.yz[k][j][i] +
                          inputs.D.b.z[k_loc] *
                                  (inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i] -
                                   inputs.E_s.xy[k + 1][j][i] - inputs.E_s.xz[k + 1][j][i]);
                } else {
                  /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                  inputs.H_s.yz[k][j][i] =
                          inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.yz[k][j][i] +
                          inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i] -
                                   inputs.E_s.xy[k + 1][j][i] - inputs.E_s.xz[k + 1][j][i]);
                }
              }
          }
          //FDTD, H_s.yz
        } else {
          //#pragma omp for
          //H_s.yz updates
          for (j = 0; j < loop_variables.J_tot_p1_bound; j++)
#pragma omp for
            for (i = 0; i < I_tot; i++) {
              for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca[n][k] = inputs.D.a.z[k_loc];
                  PSTD.cb[n][k] = inputs.D.b.z[k_loc];
                  /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                  //H_s.yz[k][j][i] = D.a.z[k_loc]*H_s.yz[k][j][i]+D.b.z[k_loc]*(E_s.xy[k][j][i] + E_s.xz[k][j][i] - E_s.xy[k+1][j][i] - E_s.xz[k+1][j][i]);
                } else {
                  PSTD.ca[n][k] = inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1];
                  PSTD.cb[n][k] = inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1];
                  /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                  //H_s.yz[k][j][i] = Dmaterial.Da.z[materials[k][j][i]-1]*H_s.yz[k][j][i]+Dmaterial.Db.z[materials[k][j][i]-1]*(E_s.xy[k][j][i] + E_s.xz[k][j][i] - E_s.xy[k+1][j][i] - E_s.xz[k+1][j][i]);
                }

                loop_variables.eh_vec[n][k][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
                loop_variables.eh_vec[n][k][1] = 0.;
              }
              k = K_tot;
              loop_variables.eh_vec[n][k][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
              loop_variables.eh_vec[n][k][1] = 0.;

              /*
    if( i==12 & j==12 ){
    for(k=0;k<K_tot;k++)
    fprintf(stdout,"%.10e ",eh_vec[n][k][0]);
    fprintf(stdout,"\n");
    }
        */

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_hz, PSTD.N_hz, inputs.H_s.yz.plan_f[n],
                               inputs.H_s.yz.plan_b[n]);

              for (k = 0; k < K_tot; k++) {
                inputs.H_s.yz[k][j][i] = PSTD.ca[n][k] * inputs.H_s.yz[k][j][i] -
                                         PSTD.cb[n][k] * loop_variables.eh_vec[n][k][0] / ((double) PSTD.N_hz);
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
              if (!inputs.materials[k][j][i]) inputs.H_s.xz[k][j][i] = 0.;
              else
                inputs.H_s.xz[k][j][i] = 0.;

#pragma omp for
        //H_s.xy update
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              k_loc = k;
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                  if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = j;
              else
                array_ind = (J_tot + 1) * k_loc + j;
              if (!inputs.materials[k][j][i])
                inputs.H_s.xy[k][j][i] =
                        inputs.D.a.y[array_ind] * inputs.H_s.xy[k][j][i] +
                        inputs.D.b.y[array_ind] *
                                (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] -
                                 inputs.E_s.zy[k][j + 1][i] - inputs.E_s.zx[k][j + 1][i]);
              else
                inputs.H_s.xy[k][j][i] =
                        inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] *
                                inputs.H_s.xy[k][j][i] +
                        inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] -
                                 inputs.E_s.zy[k][j + 1][i] - inputs.E_s.zx[k][j + 1][i]);
            }

#pragma omp for
        //H_s.yx update
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++) {
              k_loc = k;
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                  if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                      (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >=
                           (loop_variables.K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = i;
              else
                array_ind = (I_tot + 1) * k_loc + i;
              if (!inputs.materials[k][j][i])
                inputs.H_s.yx[k][j][i] =
                        inputs.D.a.x[array_ind] * inputs.H_s.yx[k][j][i] +
                        inputs.D.b.x[array_ind] *
                                (inputs.E_s.zx[k][j][i + 1] + inputs.E_s.zy[k][j][i + 1] -
                                 inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j][i]);
              else
                inputs.H_s.yx[k][j][i] =
                        inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] *
                                inputs.H_s.yx[k][j][i] +
                        inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                (inputs.E_s.zx[k][j][i + 1] + inputs.E_s.zy[k][j][i + 1] -
                                 inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j][i]);
            }

#pragma omp for
        for (k = 0; k <= K_tot; k++) {
          for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++)
              if (!inputs.materials[k][j][i]) inputs.H_s.yz[k][j][i] = 0.;
              else
                inputs.H_s.yz[k][j][i] = 0.;
        }
      }

      if (inputs.params.dimension == THREE ||
          inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
//FDTD, H_s.zy
#pragma omp for
          //H_s.zy update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.zy[k][j][i] =
                          inputs.D.a.y[array_ind] * inputs.H_s.zy[k][j][i] +
                          inputs.D.b.y[array_ind] *
                                  (inputs.E_s.xy[k][j + 1][i] + inputs.E_s.xz[k][j + 1][i] -
                                   inputs.E_s.xy[k][j][i] - inputs.E_s.xz[k][j][i]);
                else
                  inputs.H_s.zy[k][j][i] =
                          inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.zy[k][j][i] +
                          inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.xy[k][j + 1][i] + inputs.E_s.xz[k][j + 1][i] -
                                   inputs.E_s.xy[k][j][i] - inputs.E_s.xz[k][j][i]);
              }
          //FDTD, H_s.zy
        } else {
#pragma omp for
          //H_s.zy update
          for (k = 0; k < (K_tot + 1); k++)
            for (i = 0; i < I_tot; i++) {
              for (j = 0; j < J_tot; j++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i]) {
                  PSTD.ca[n][j] = inputs.D.a.y[array_ind];
                  PSTD.cb[n][j] = inputs.D.b.y[array_ind];
                  //	      H_s.zy[k][j][i] = D.a.y[array_ind]*H_s.zy[k][j][i]+D.b.y[array_ind]*(E_s.xy[k][j+1][i] + E_s.xz[k][j+1][i] - E_s.xy[k][j][i] - E_s.xz[k][j][i]);
                } else {
                  PSTD.ca[n][j] = inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1];
                  PSTD.cb[n][j] = inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1];
                  //	      H_s.zy[k][j][i] = Dmaterial.Da.y[materials[k][j][i]-1]*H_s.zy[k][j][i]+Dmaterial.Db.y[materials[k][j][i]-1]*(E_s.xy[k][j+1][i] + E_s.xz[k][j+1][i] - E_s.xy[k][j][i] - E_s.xz[k][j][i]);
                }

                loop_variables.eh_vec[n][j][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
                loop_variables.eh_vec[n][j][1] = 0.;
              }
              j = J_tot;
              loop_variables.eh_vec[n][j][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
              loop_variables.eh_vec[n][j][1] = 0.;

              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_hy, PSTD.N_hy, inputs.H_s.zy.plan_f[n],
                               inputs.H_s.zy.plan_b[n]);

              for (j = 0; j < J_tot; j++) {
                inputs.H_s.zy[k][j][i] = PSTD.ca[n][j] * inputs.H_s.zy[k][j][i] +
                                         PSTD.cb[n][j] * loop_variables.eh_vec[n][j][0] / ((double) PSTD.N_hy);
              }
            }
          //PSTD, H_s.zy
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)


        if (solver_method == SolverMethod::FiniteDifference) {
//FDTD, H_s.zx
#pragma omp for
          //H_s.zx update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_tot_bound; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.zx[k][j][i] =
                          inputs.D.a.x[array_ind] * inputs.H_s.zx[k][j][i] +
                          inputs.D.b.x[array_ind] *
                                  (inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i] -
                                   inputs.E_s.yx[k][j][i + 1] - inputs.E_s.yz[k][j][i + 1]);
                else
                  inputs.H_s.zx[k][j][i] =
                          inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] *
                                  inputs.H_s.zx[k][j][i] +
                          inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                  (inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i] -
                                   inputs.E_s.yx[k][j][i + 1] - inputs.E_s.yz[k][j][i + 1]);
              }
          //FDTD, H_s.zx
        } else {
#pragma omp for
          //H_s.zx update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < loop_variables.J_tot_bound; j++) {
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + loop_variables.K)) {
                    if ((k - inputs.structure[i][1]) < (loop_variables.K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >=
                             (loop_variables.K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + loop_variables.K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i]) {
                  //		H_s.zx[k][j][i] = D.a.x[array_ind]*H_s.zx[k][j][i]+D.b.x[array_ind]*(E_s.yx[k][j][i] + E_s.yz[k][j][i] - E_s.yx[k][j][i+1] - E_s.yz[k][j][i+1]);
                  PSTD.ca[n][i] = inputs.D.a.x[array_ind];
                  PSTD.cb[n][i] = inputs.D.b.x[array_ind];
                } else {
                  //		H_s.zx[k][j][i] = Dmaterial.Da.x[materials[k][j][i]-1]*H_s.zx[k][j][i]+Dmaterial.Db.x[materials[k][j][i]-1]*(E_s.yx[k][j][i] + E_s.yz[k][j][i] - E_s.yx[k][j][i+1] - E_s.yz[k][j][i+1]);
                  PSTD.ca[n][i] = inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1];
                  PSTD.cb[n][i] = inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1];
                }

                loop_variables.eh_vec[n][i][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
                loop_variables.eh_vec[n][i][1] = 0.;
              }
              i = I_tot;
              loop_variables.eh_vec[n][i][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
              loop_variables.eh_vec[n][i][1] = 0.;


              first_derivative(loop_variables.eh_vec[n], loop_variables.eh_vec[n], PSTD.dk_hx, PSTD.N_hx, inputs.H_s.zx.plan_f[n],
                               inputs.H_s.zx.plan_b[n]);

              for (i = 0; i < I_tot; i++) {
                inputs.H_s.zx[k][j][i] = PSTD.ca[n][i] * inputs.H_s.zx[k][j][i] -
                                         PSTD.cb[n][i] * loop_variables.eh_vec[n][i][0] / ((double) PSTD.N_hx);
              }
            }
          //PSTD, H_s.zx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
      }  //(params.dimension==THREE || params.dimension==TE)
    }    //end parallel
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    //fprintf(stderr,"Pos 11b:\n");
    //update terms for self consistency across scattered/total interface - E updates
    if (inputs.params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase =
              exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_E, 2. * DCPI));
      double commonAmplitude = linear_ramp(time_E);
      for (k = (inputs.K0.index); k <= (inputs.K1.index); k++)
        for (j = (inputs.J0.index); j <= (inputs.J1.index); j++) {
          if (inputs.I0.apply) {//Perform across I0

            if (!inputs.params.is_multilayer) array_ind = inputs.I0.index - 1;
            else
              array_ind = (I_tot + 1) * k + inputs.I0.index - 1;

            if (j < (inputs.J1.index))
              inputs.H_s.zx[k][j][(inputs.I0.index) - 1] =
                      inputs.H_s.zx[k][j][(inputs.I0.index) - 1] +
                      inputs.D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][0] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][0]));
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.yx[k][j][(inputs.I0.index) - 1] =
                      inputs.H_s.yx[k][j][(inputs.I0.index) - 1] -
                      inputs.D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][1] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][1]));
          }
          if (inputs.I1.apply) {//Perform across I1

            if (!inputs.params.is_multilayer) array_ind = inputs.I1.index;
            else
              array_ind = (I_tot + 1) * k + inputs.I1.index;

            if (j < (inputs.J1.index))
              inputs.H_s.zx[k][j][(inputs.I1.index)] =
                      inputs.H_s.zx[k][j][(inputs.I1.index)] -
                      inputs.D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][4] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][4]));
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.yx[k][j][(inputs.I1.index)] =
                      inputs.H_s.yx[k][j][(inputs.I1.index)] +
                      inputs.D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource
                                            .real[k - (inputs.K0.index)][j - (inputs.J0.index)][5] +
                                    IMAGINARY_UNIT *
                                            inputs.Isource.imag[k - (inputs.K0.index)]
                                                               [j - (inputs.J0.index)][5]));
          }
        }

      for (k = (inputs.K0.index); k <= (inputs.K1.index); k++)
        for (i = (inputs.I0.index); i <= (inputs.I1.index); i++) {
          if (inputs.J0.apply) {//Perform across J0

            if (!inputs.params.is_multilayer) array_ind = inputs.J0.index;
            else
              array_ind = (J_tot + 1) * k + inputs.J0.index;

            if (i < (inputs.I1.index))
              inputs.H_s.zy[k][(inputs.J0.index) - 1][i] =
                      inputs.H_s.zy[k][(inputs.J0.index) - 1][i] -
                      inputs.D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][0] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][0]));

            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.xy[k][(inputs.J0.index) - 1][i] =
                      inputs.H_s.xy[k][(inputs.J0.index) - 1][i] +
                      inputs.D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][1] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][1]));
          }
          if (inputs.J1.apply) {//Perform across J1

            if (!inputs.params.is_multilayer) array_ind = inputs.J1.index;
            else
              array_ind = (J_tot + 1) * k + inputs.J1.index;

            if (i < (inputs.I1.index))
              inputs.H_s.zy[k][(inputs.J1.index)][i] =
                      inputs.H_s.zy[k][(inputs.J1.index)][i] +
                      inputs.D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][4] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][4]));
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.xy[k][(inputs.J1.index)][i] =
                      inputs.H_s.xy[k][(inputs.J1.index)][i] -
                      inputs.D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource
                                            .real[k - (inputs.K0.index)][i - (inputs.I0.index)][5] +
                                    IMAGINARY_UNIT *
                                            inputs.Jsource.imag[k - (inputs.K0.index)]
                                                               [i - (inputs.I0.index)][5]));
          }
        }

      for (j = (inputs.J0.index); j <= (inputs.J1.index); j++)
        for (i = (inputs.I0.index); i <= (inputs.I1.index); i++) {
          if (inputs.K0.apply) {//Perform across K0
            if (i < (inputs.I1.index))
              inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                      inputs.D.b.z[(inputs.K0.index) - 1] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][0] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][0]));
            if (j < (inputs.J1.index))
              inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.xz[(inputs.K0.index) - 1][j][i] -
                      inputs.D.b.z[(inputs.K0.index) - 1] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][1] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][1]));
          }
          if (inputs.K1.apply) {//Perform across K1
            if (i < (inputs.I1.index))
              inputs.H_s.yz[(inputs.K1.index)][j][i] =
                      inputs.H_s.yz[(inputs.K1.index)][j][i] -
                      inputs.D.b.z[(inputs.K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][4] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][4]));
            if (j < (inputs.J1.index))
              inputs.H_s.xz[(inputs.K1.index)][j][i] =
                      inputs.H_s.xz[(inputs.K1.index)][j][i] +
                      inputs.D.b.z[(inputs.K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource
                                            .real[j - (inputs.J0.index)][i - (inputs.I0.index)][5] +
                                    IMAGINARY_UNIT *
                                            inputs.Ksource.imag[j - (inputs.J0.index)]
                                                               [i - (inputs.I0.index)][5]));
          }
        }
      outputs.E.ft = real(commonAmplitude * commonPhase);
    } else if (inputs.params.source_mode == 1) {//pulsed
      //fprintf(stderr,"Pos 11c\n");
      if (J_tot == 0) {
        //fprintf(stderr,"Pos 11d\n");
        j = 0;
        for (i = 0; i < (I_tot + 1); i++) {
          inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                  inputs.H_s.xz[(inputs.K0.index) - 1][j][i] -
                  inputs.D.b.z[(inputs.K0.index) - 1] *
                          real((inputs.Ksource.real[0][i - (inputs.I0.index)][1] +
                                IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][1]) *
                               (-1. * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT *
                                   fmod(inputs.params.omega_an * (time_E - inputs.params.to_l),
                                        2 * DCPI))) *
                          exp(-1. * DCPI *
                              pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
          //broadband source term
          if (inputs.params.eyi_present)
            inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.xz[(inputs.K0.index) - 1][j][i] -
                    inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.y[tind][j][i];
        }
        //fprintf(stderr,"Pos 11e\n");
        for (i = 0; i < I_tot; i++) {
          inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                  inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                  inputs.D.b.z[(inputs.K0.index) - 1] *
                          real((inputs.Ksource.real[0][i - (inputs.I0.index)][0] +
                                IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][0]) *
                               (-1. * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT *
                                   fmod(inputs.params.omega_an * (time_E - inputs.params.to_l),
                                        2 * DCPI))) *
                          exp(-1. * DCPI *
                              pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
          //broadband source term
          if (inputs.params.exi_present)
            inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                    inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.x[tind][j][i];
          //if(i==511)
          //  fprintf(stdout,"%e\n",D.b.z[((int)K0[0])-1]*exi[tind][j][i]);
        }
        //fprintf(stderr,"Pos 11f\n");
      } else {
        //fprintf(stderr,"Pos 11g\n");
        for (j = 0; j < J_tot; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.xz[(inputs.K0.index) - 1][j][i] -
                    inputs.D.b.z[(inputs.K0.index) - 1] *
                            real((inputs.Ksource
                                          .real[j - (inputs.J0.index)][i - (inputs.I0.index)][1] +
                                  IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                      [i - (inputs.I0.index)][1]) *
                                 (-1. * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT *
                                     fmod(inputs.params.omega_an * (time_E - inputs.params.to_l),
                                          2 * DCPI))) *
                            exp(-1. * DCPI *
                                pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
            //broadband source term
            if (inputs.params.eyi_present)
              inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.xz[(inputs.K0.index) - 1][j][i] -
                      inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.y[tind][j][i];
          }
        //fprintf(stderr,"Pos 11h\n");
        for (j = 0; j < (J_tot + 1); j++)
          for (i = 0; i < I_tot; i++) {
            inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                    inputs.D.b.z[(inputs.K0.index) - 1] *
                            real((inputs.Ksource
                                          .real[j - (inputs.J0.index)][i - (inputs.I0.index)][0] +
                                  IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)]
                                                                      [i - (inputs.I0.index)][0]) *
                                 (-1. * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT *
                                     fmod(inputs.params.omega_an * (time_E - inputs.params.to_l),
                                          2 * DCPI))) *
                            exp(-1. * DCPI *
                                pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
            //broadband source term
            if (inputs.params.exi_present)
              inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                      inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.x[tind][j][i];
          }
        //fprintf(stderr,"Pos 11i\n");
      }
      outputs.E.ft =
              real((-1. * IMAGINARY_UNIT) *
                   exp(-IMAGINARY_UNIT *
                       fmod(inputs.params.omega_an * (time_E - inputs.params.to_l), 2 * DCPI))) *
              exp(-1. * DCPI * pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
      //fprintf(stderr,"Pos 11j\n");
    }
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    if (inputs.params.exphasorssurface || inputs.params.exphasorsvolume ||
        inputs.params.exdetintegral || outputs.vertex_phasors.there_are_vertices_to_extract_at()) {
      if (inputs.params.source_mode == SourceMode::steadystate) {
        /*
	  Each time a new acquisition period of harmonic illumination begins, all complex amplitudes
	  (volume, surface etc.) are set back to 0 since the discrete Fourier transforms used to acquire
	  these complex amplitudes starts again. In particular, the returned complex amplitudes will have
	  been acquired during a single acquisition period of harmonic illumination. Note that, as
	  explained above, the acquisition period is actually three periods of the harmonic waves
	  fundamental period. The complex amplitudes are reset to 0 using calls such as:

    However, the normalisation factors are reset to 0 here.
	 */

        if ((tind % inputs.Nsteps) == 0) {
          outputs.E.angular_norm = 0.0;
          outputs.H.angular_norm = 0.0;

          for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
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

        outputs.E.add_to_angular_norm(tind, inputs.Nsteps, inputs.params);
        outputs.H.add_to_angular_norm(tind, inputs.Nsteps, inputs.params);

        for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
          extract_phasor_norms(ifx, tind, inputs.Nsteps);
        }
      } else {
        if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {

          outputs.E.add_to_angular_norm(tind, inputs.params.Npe, inputs.params);
          outputs.H.add_to_angular_norm(tind, inputs.params.Npe, inputs.params);

          for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
            extract_phasor_norms(ifx, tind, inputs.params.Npe);
          }
        }
      }
    }
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    if ((((double) time(NULL)) - t0) > 1) {

      maxfield = max(inputs.E_s.largest_field_value(), inputs.H_s.largest_field_value());

      spdlog::info("Iterating: tind = {0:d}, maxfield = {1:e}", tind, maxfield);
      t0 = double(time(NULL));
    }
    //fprintf(stderr,"Post-iter 3\n");
    if ((inputs.params.source_mode == SourceMode::steadystate) &&
        (tind == (inputs.params.Nt - 1)) && (inputs.params.run_mode == RunMode::complete) &&
        inputs.params.exphasorsvolume) {
      fprintf(stdout, "Iteration limit reached, setting output fields to last complete DFT\n");
      outputs.E.set_values_from(loop_variables.E_copy);
    }
    //fprintf(stderr,"Post-iter 4\n");
    fflush(stdout);
    //fprintf(stderr,"Post-iter 5\n");
    //fprintf(stderr,"%s %d %d\n",tdfdirstr, strcmp(tdfdirstr,""),are_equal(tdfdirstr,""));
    if (inputs.params.has_tdfdir && (tind % inputs.params.Np) == 0) {
      fprintf(stderr, "Saving field\n");
      inputs.ex_td_field_exporter.export_field(inputs.E_s, inputs.skip_tdf, tind);
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
    timers.end_timer(TimersTrackingLoop::MAIN);
    spdlog::info("Time elapsed in main loop (s): {0:e}",
                 timers.time_ellapsed_by(TimersTrackingLoop::MAIN));
  }
}