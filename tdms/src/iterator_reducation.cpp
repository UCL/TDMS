#include "iterator.h"

#include <complex>
#include <algorithm>

#include <omp.h>
#include <spdlog/spdlog.h>

#include "mesh_base.h"
#include "numerical_derivative.h"
#include "array_init.h"
#include "fieldsample.h"
#include "globals.h"
#include "interface.h"
#include "iterator_intermediate_matlab_variables.h"
#include "iterator_loop_variables.h"
#include "iterator_main_loop.h"
#include "matlabio.h"
#include "shapes.h"
#include "source.h"
#include "surface_phasors.h"
#include "vertex_phasors.h"
#include "utils.h"

using namespace std;
using namespace tdms_math_constants;
using namespace tdms_phys_constants;

//whether of not to time execution
#define TIME_EXEC false
//time main loop
#define TIME_MAIN_LOOP true
//threshold used to terminate the steady state iterations
#define TOL 1e-6
//parameter to control the PreferredInterpolationMethodswith of the ramp when introducing the waveform in steady state mode
#define ramp_width 4.

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
void execute_simulation(int nlhs, mxArray *plhs[], int nrhs, InputMatrices in_matrices,
                        SolverMethod solver_method,
                        PreferredInterpolationMethods preferred_interpolation_methods) {

  if (solver_method == SolverMethod::FiniteDifference) {
    spdlog::info("Using finite-difference method (FDTD)");
  } else {
    spdlog::info("Using pseudospectral method (PSTD)");
  }
  if (preferred_interpolation_methods == PreferredInterpolationMethods::BandLimited) {
    spdlog::info("Using band-limited interpolation where possible");
  } else {
    spdlog::info("Restricting to cubic interpolation");
  }

  Iterator_LoopVariables main_loop_variables(in_matrices, solver_method, preferred_interpolation_methods);

  double maxfield = 0.;

  IteratorMainLoop iteration_controller;

  // mwSize *dims;
  // dims = (mwSize *) malloc(3 * sizeof(mwSize));
  // mwSize *label_dims;
  // label_dims = (mwSize *) malloc(2 * sizeof(mwSize));
  mxArray *mx_surface_facets;//< surface_facets RECYCLED after the main loop for the outputs, but isn't actually needed at this scope for the setup and main loop

  Iterator_IntermediateMATLABVariables iMVars;

  spdlog::info("Using {} OMP threads", omp_get_max_threads());

  if (nrhs != 49) { throw runtime_error("Expected 49 inputs. Had " + to_string(nrhs)); }

  if (nlhs != 31) { throw runtime_error("31 outputs required. Had " + to_string(nlhs)); }

  // I_LoopVariable beings handling from here
  main_loop_variables.link_fields_and_labels_to_output_pointers(plhs);
  //plhs[13] -> plhs[15] are the interpolated electric field values
  //plhs[16] -> plhs[18] are the interpolated magnetic field values

  main_loop_variables.link_id_to_output_pointers(plhs, iMVars);

  /*Now set up the phasor arrays for storing the fdtd version of the input fields,
    these will be used in a boot strapping procedure. Calculated over a complete
    xy-plane. */

  // TO HERE WITH REDUCTION

  ndims = 2;
  dims[0] = I_tot;
  dims[1] = J_tot + 1;
  plhs[6] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                 mxCOMPLEX);//x electric field source phasor - boot strapping
  iwave_lEx_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) plhs[6]), dims[0], dims[1]);
  iwave_lEx_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) plhs[6]), dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lEx_Rbs, dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lEx_Ibs, dims[0], dims[1]);

  dims[0] = I_tot + 1;
  dims[1] = J_tot;
  plhs[7] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                 mxCOMPLEX);//y electric field source phasor - boot strapping
  iwave_lEy_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) plhs[7]), dims[0], dims[1]);
  iwave_lEy_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) plhs[7]), dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lEy_Rbs, dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lEy_Ibs, dims[0], dims[1]);

  dims[0] = I_tot + 1;
  dims[1] = J_tot;
  plhs[8] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                 mxCOMPLEX);//x magnetic field source phasor - boot strapping

  iwave_lHx_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) plhs[8]), dims[0], dims[1]);
  iwave_lHx_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) plhs[8]), dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lHx_Rbs, dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lHx_Ibs, dims[0], dims[1]);

  dims[0] = I_tot;
  dims[1] = J_tot + 1;
  plhs[9] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                 mxCOMPLEX);//y magnetic field source phasor - boot strapping
  iwave_lHy_Rbs = cast_matlab_2D_array(mxGetPr((mxArray *) plhs[9]), dims[0], dims[1]);
  iwave_lHy_Ibs = cast_matlab_2D_array(mxGetPi((mxArray *) plhs[9]), dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lHy_Rbs, dims[0], dims[1]);
  initialiseDouble2DArray(iwave_lHy_Ibs, dims[0], dims[1]);

  /*start dispersive*/

  //work out if we have any disperive materials
  bool is_disp = is_dispersive(materials, gamma, params.dt, I_tot, J_tot, K_tot);
  //work out if we have conductive background: background is conductive if at least one entry exceeds 1e-15
  bool is_conductive = !(rho_cond.all_elements_less_than(1e-15, I_tot + 1, J_tot + 1, K_tot + 1));
  // work out if we have a dispersive background
  if (params.is_disp_ml) params.is_disp_ml = matched_layer.is_dispersive(K_tot);
  //  fprintf(stderr,"is_disp:%d, is_conductive%d, params.is_disp_ml: %d\n",is_disp,is_conductive,params.is_disp_ml);
  //if we have dispersive materials we need to create additional field variables
  auto E_nm1 = ElectricSplitField(I_tot, J_tot, K_tot);
  auto J_nm1 = CurrentDensitySplitField(I_tot, J_tot, K_tot);

  if (is_disp || params.is_disp_ml) {
    E_nm1.allocate_and_zero();
    J_nm1.allocate_and_zero();
    J_s.allocate_and_zero();
  }
  //fprintf(stderr,"Pre 14\n");
  auto J_c = CurrentDensitySplitField(I_tot, J_tot, K_tot);
  if (is_conductive) { J_c.allocate_and_zero(); }
  /*end dispersive*/

  plhs[27] = fieldsample.mx;

  vertex_phasors.setup_camplitude_arrays(f_ex_vec.size());
  plhs[28] = vertex_phasors.get_mx_camplitudes();

  /*end of setup the output array for the sampled field*/


  /*set up the parameters for the phasor convergence procedure*/
  /*First we set dt so that an integer number of time periods fits within a sinusoidal period
   */
  double Nsteps_tmp = 0.0;
  double dt_old;
  if (params.source_mode == SourceMode::steadystate) {
    dt_old = params.dt;
    Nsteps_tmp = ceil(2. * DCPI / params.omega_an / params.dt * 3);
    params.dt = 2. * DCPI / params.omega_an * 3 / Nsteps_tmp;
  }

  //fprintf(stderr,"Pre 16\n");
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete)
    fprintf(stderr, "Changing dt from %.10e to %.10e\n", dt_old, params.dt);
  iteration_controller.Nsteps = (int) lround(Nsteps_tmp);
  //fprintf(stderr,"Pre 17\n");
  //iteration_controller.Nsteps = (int)(floor(3*2.*DCPI/(params.omega_an*params.dt)) + 1.);//the number of time steps in a sinusoidal period
  //fprintf(stderr,"Pre 18\n");
  /*params.Nt should be an integer number of iteration_controller.Nsteps in the case of steady-state operation*/
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete)
    if (params.Nt / iteration_controller.Nsteps * iteration_controller.Nsteps != params.Nt) {
      fprintf(stderr, "Changing the value of Nt from %d to", params.Nt);
      params.Nt = params.Nt / iteration_controller.Nsteps * iteration_controller.Nsteps;
      fprintf(stderr, " %d for correct phasor extraction\n", params.Nt);
    }
  //fprintf(stderr,"Pre 19\n");

  if ((params.run_mode == RunMode::complete) && (params.source_mode == SourceMode::steadystate)) printf("iteration_controller.Nsteps: %d \n", iteration_controller.Nsteps);

  /*An optimization step in the 2D (J_tot==0) case, try to work out if we have either
    of TE or TM, ie, not both*/
  int ksource_nz[4];
  for (int icomp = 0; icomp < 4; icomp++) ksource_nz[icomp] = 0;
  //fprintf(stderr,"Pre 20\n");
  if (J_tot == 0) {


    for (int icomp = 0; icomp < 4; icomp++)
      for (i = 0; i < (I_tot + 1); i++) {
        ksource_nz[icomp] = ksource_nz[icomp] ||
                            (fabs(Ksource.imag[0][i - (I0.index)][icomp]) > 1.0e-15) ||
                            (fabs(Ksource.real[0][i - (I0.index)][icomp]) > 1.0e-15);
      }

    //for (int icomp=0;icomp<4;icomp++)
    //	fprintf(stderr,"ksource_nz[%d] = %d\n",icomp,ksource_nz[icomp]);
  }
  //fprintf(stderr,"Pre 21\n");
  /*
    In the J_tot==0 2D version, the 'TE' case involves components Ey, Hx and Hz
    'TM' case involves components Ex, Ez and Hy

    Ey and Hx receive an input from the source condition only if ksource_nz[2] or ksource_nz[1]
    are non-zero.

    Ex and Hy receive an input from the source condition only if ksource_nz[3] or ksource_nz[0]
    are non-zero.

    The idea is to use an alternative upper limit to the loop over j when we have J_tot==0. We currently have the followingloops
    in place for each of the component updates.

    From the analysis below, we see that in all cases, the TE update has the following loop on j:
    for(j=0;j<max(J_tot,1);j++)
    whilst the TM case has:
    for(j=0;j<(J_tot+1);j++)

    So we can create variables
    J_tot_p1_bound
    and
    J_tot_bound

    which would take the following values
    3D:
    J_tot_p1_bound = J_tot + 1;
    J_tot_bound = J_tot;

    2D:
    TE:
    J_tot_bound = 1;
    Not TE:
    J_tot_bound = 0;

    TM:
    J_tot_p1_bound = 1;
    Not TM:
    J_tot_p1_bound = 0;

    Exy: Not involved in 2D
    for(k=0;k<(K_tot+1);k++)
    for(j=1;j<J_tot;j++)
    for(i=0;i<I_tot;i++){

    Exz: TM
    for(k=1;k<K_tot;k++)
    for(j=0;j<(J_tot+1);j++)
    for(i=0;i<I_tot;i++){

    Eyx: TE
    for(k=0;k<(K_tot+1);k++)
    for(j=0;j<max(J_tot,1);j++)
    for(i=1;i<I_tot;i++){

    Eyz: TE
    for(k=1;k<K_tot;k++)
    for(j=0;j<max(J_tot,1);j++)
    for(i=0;i<(I_tot+1);i++){

    Ezx: TM
    for(k=0;k<K_tot;k++)
    for(j=0;j<(J_tot+1);j++)
    for(i=1;i<I_tot;i++){

    Ezy: Not involved in 2D
    for(k=0;k<K_tot;k++)
    for(j=1;j<J_tot;j++)
    for(i=0;i<(I_tot+1);i++){

    Hxy:  Not involved in 2D
    for(k=0;k<K_tot;k++)
    for(j=0;j<J_tot;j++)
    for(i=0;i<(I_tot+1);i++)

    Hxz: TE
    for(k=0;k<K_tot;k++)
    for(j=0;j<max(J_tot,1);j++)
    for(i=0;i<(I_tot+1);i++){

    Hyx: TM
    for(k=0;k<K_tot;k++)
    for(j=0;j<(J_tot+1);j++)
    for(i=0;i<I_tot;i++){

    Hyz: TM
    for(k=0;k<K_tot;k++){
    for(j=0;j<(J_tot+1);j++)
    for(i=0;i<I_tot;i++){

    Hzx: TE
    for(k=0;k<(K_tot+1);k++)
    for(j=0;j<max(J_tot,1);j++)
    for(i=0;i<I_tot;i++){

    Hzy: Not involved in 2D
    for(k=0;k<(K_tot+1);k++)
    for(j=0;j<J_tot;j++)
    for(i=0;i<I_tot;i++){
  */
  //implement the above
  int J_tot_bound = J_tot;
  int J_tot_p1_bound = J_tot + 1;
  if (J_tot == 0) {
    //TE case
    if (ksource_nz[2] || ksource_nz[1] || params.eyi_present) J_tot_bound = 1;
    else
      J_tot_bound = 0;

    //TM case
    if (ksource_nz[3] || ksource_nz[0] || params.exi_present) J_tot_p1_bound = 1;
    else
      J_tot_p1_bound = 0;
  }
  //fprintf(stderr,"Pre 22a\n");
  //fprintf(stderr,"Pre 23\n");

  /*Start of FDTD iteration*/
  //open a file for logging the times

  spdlog::debug("Starting main loop");

  if (TIME_MAIN_LOOP) { iteration_controller.start_timer(IterationTimers::MAIN); }

  for (iteration_controller.tind = params.start_tind; iteration_controller.tind < params.Nt; iteration_controller.tind++) {
    // Update the "time" the fields are currently at
    iteration_controller.update_field_times_to_current_iteration(params.dt);

    //Extract phasors
    iteration_controller.start_timer(IterationTimers::INTERNAL);
    if ((iteration_controller.dft_counter == iteration_controller.Nsteps) && (params.run_mode == RunMode::complete)
        && (params.source_mode == SourceMode::steadystate) && params.exphasorsvolume) {

      iteration_controller.dft_counter = 0;

      double tol = E.normalised_difference(E_copy);
      if (tol < TOL) break; //required accuracy obtained

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

    if ((params.source_mode == SourceMode::steadystate) && (params.run_mode == RunMode::complete) && params.exphasorsvolume) {

      E.set_phasors(E_s, iteration_controller.dft_counter - 1, params.omega_an, params.dt, iteration_controller.Nsteps);
      H.set_phasors(H_s, iteration_controller.dft_counter, params.omega_an, params.dt, iteration_controller.Nsteps);

      if (params.exphasorssurface) {
        if (params.intphasorssurface) {
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            surface_phasors.extractPhasorsSurface(ifx, E_s, H_s,
                                                  iteration_controller.dft_counter, f_ex_vec[ifx] * 2 * DCPI, iteration_controller.Nsteps,
                                                  params);
          }
          iteration_controller.dft_counter++;
        } else {
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            // do not interpolate when extracting
            surface_phasors.extractPhasorsSurface(ifx, E_s, H_s, iteration_controller.dft_counter,
                    f_ex_vec[ifx] * 2 * DCPI, iteration_controller.Nsteps, params, false);
          }
          iteration_controller.dft_counter++;
        }
      }

    } else if ((params.source_mode == SourceMode::pulsed) && (params.run_mode == RunMode::complete) && params.exphasorsvolume) {
      if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }

      if ((iteration_controller.tind - params.start_tind) % params.Np == 0) {
        E.set_phasors(E_s, iteration_controller.tind - 1, params.omega_an, params.dt, params.Npe);
        H.set_phasors(H_s, iteration_controller.tind, params.omega_an, params.dt, params.Npe);
      }
      if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }
      //fprintf(stderr,"Pos 01b:\n");
    }

    /*extract fieldsample*/
    if (fieldsample.all_vectors_are_non_empty()) {
      fieldsample.extract(E_s, params.pml, params.Nt);
    }
    /*end extract fieldsample*/

    //fprintf(stderr,"Pos 02:\n");
    if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete && params.exphasorssurface) {
      if ((iteration_controller.tind - params.start_tind) % params.Np == 0) {
        if (params.intphasorssurface)
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
              surface_phasors.extractPhasorsSurface(ifx, E_s, H_s,
                                                    iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI, params.Npe,
                                                    params);
          }
        else
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            // do not interpolate when extracting
              surface_phasors.extractPhasorsSurface(ifx, E_s, H_s, iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI,
                      params.Npe, params, false);
          }
      }
    }

    if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete &&
        (vertex_phasors.there_are_vertices_to_extract_at()) &&
        ((iteration_controller.tind - params.start_tind) % params.Np == 0)) {
      //     fprintf(stderr,"loc 01 (%d,%d,%d)\n",iteration_controller.tind,params.start_tind,params.Np);
          //fprintf(stderr,"loc 03\n");
          //	  fprintf(stderr,"EPV 01\n");
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            vertex_phasors.extractPhasorsVertices(ifx, E_s, H_s, iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI, params);
          }
    }

    //fprintf(stderr,"Pos 02a:\n");
    if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete && params.exdetintegral) {
      if ((iteration_controller.tind - params.start_tind) % params.Np == 0) {
        spdlog::debug("Setting Ex_t, Ey_t");

        //First need to sum up the Ex and Ey values on a plane ready for FFT, remember that Ex_t and
        // Ey_t are in row-major format whilst Exy etc. are in column major format
        for (j = params.pml.Dyl; j < (J_tot - params.pml.Dyu); j++)
          for (i = params.pml.Dxl; i < (I_tot - params.pml.Dxu); i++) {
            int m = j - params.pml.Dyl + (i - params.pml.Dxl) * (J_tot - params.pml.Dyu - params.pml.Dyl);
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
                      kprop = exp(IMAGINARY_UNIT * (-params.air_interface + params.z_obs) * 2. * DCPI / lambda_an_t * refind *
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
              phaseTermE = fmod(f_ex_vec[ifx] * 2. * DCPI * ((double) iteration_controller.tind) * params.dt, 2 * DCPI);
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
        extractPhasorsPlane(iwave_lEx_Rbs, iwave_lEx_Ibs, iwave_lEy_Rbs, iwave_lEy_Ibs,
                            iwave_lHx_Rbs, iwave_lHx_Ibs, iwave_lHy_Rbs, iwave_lHy_Ibs, E_s, H_s,
                            I_tot, J_tot, K0.index + 1, iteration_controller.tind, params.omega_an, params.dt,
                            params.Nt);//extract the phasors just above the line
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
    if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }
    //fprintf(stderr,"Dimension = %d\n",params.dimension);
    /*
      for(k=0;k<(K_tot+1);k++)
      fprintf(stdout,"%e ",Exy[k][13][13]+Exz[k][13][13]);
      fprintf(stdout,"\n");
    */
    (void)n;  // n is unused in FD derivatives – this silences the compiler warning

#pragma omp parallel default(shared) private(i, j, k, n, rho, k_loc, array_ind, Ca, Cb, Cc, alpha_l,\
                                             beta_l, gamma_l, kappa_l, sigma_l, Enp1,              \
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
                  //fprintf(stdout,"(%d,%d,%d,%d)\n",i,j,k,iteration_controller.tind);
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
                  //fprintf(stdout,"(%d,%d,%d,%d)\n",i,j,k,iteration_controller.tind);
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
        } // if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
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
      if(iteration_controller.tind==0)
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
      if(iteration_controller.tind==0)
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
    if (iteration_controller.tind==1 & i==25 & j==25){
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
    if (iteration_controller.tind==1 & i==25 & j==25){
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
        } // if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)

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
      }//(params.dimension==THREE || params.dimension==TE)
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
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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
      }//(params.dimension==THREE || params.dimension==TE)
      else {
#pragma omp for
        for (k = 0; k <= K_tot; k++)
          for (j = 1; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              rho = 0.;
              k_loc = k;
              if (params.is_structure)
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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
    if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }
    /********************/

    //update terms for self consistency across scattered/total interface - E updates##
    if (params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase = exp(-IMAGINARY_UNIT * fmod(params.omega_an * iteration_controller.time_H, 2. * DCPI));
      double commonAmplitude = linearRamp(iteration_controller.time_H, 1. / (params.omega_an / (2 * DCPI)), ramp_width);
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
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
              if (is_conductive)
                J_c.zx[k][j][I0.index] +=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                              IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
              if (params.is_disp_ml)
                J_s.zx[k][j][I0.index] +=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
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
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
              if (is_conductive)
                J_c.yx[k][j][I0.index] -=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                              IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
              if (params.is_disp_ml)
                J_s.yx[k][j][I0.index] -=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
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
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
              if (is_conductive)
                J_c.zx[k][j][I1.index] -=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                              IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
              if (params.is_disp_ml)
                J_s.zx[k][j][I1.index] -=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
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
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
              if (is_conductive)
                J_c.yx[k][j][I1.index] +=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                              IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
              if (params.is_disp_ml)
                J_s.yx[k][j][I1.index] +=
                        matched_layer.kappa.x[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
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
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
              if (is_conductive)
                J_c.zy[k][(J0.index)][i] -=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                              IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
              if (params.is_disp_ml)
                J_s.zy[k][(J0.index)][i] -=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
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
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
              if (is_conductive)
                J_c.xy[k][(J0.index)][i] +=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                              IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
              if (params.is_disp_ml)
                J_s.xy[k][(J0.index)][i] +=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
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
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
              if (is_conductive)
                J_c.zy[k][(J1.index)][i] +=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                              IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
              if (params.is_disp_ml)
                J_s.zy[k][(J1.index)][i] -=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
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
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
              if (is_conductive)
                J_c.xy[k][(J1.index)][i] -=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                              IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
              if (params.is_disp_ml)
                J_s.xy[k][(J1.index)][i] +=
                        matched_layer.kappa.y[array_ind] * matched_layer.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
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
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
              if (is_conductive)
                J_c.yz[(K0.index)][j][i] +=
                        rho_cond.z[(K0.index)] * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                              IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
              if (params.is_disp_ml)
                J_s.yz[(K0.index)][j][i] -=
                        matched_layer.kappa.z[(K0.index)] * matched_layer.gamma[k] / (2. * params.dt) * C.b.z[K0.index] *
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
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
              if (is_conductive)
                J_c.xz[(K0.index)][j][i] -=
                        rho_cond.z[(K0.index)] * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                              IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
              if (params.is_disp_ml)
                J_s.xz[(K0.index)][j][i] +=
                        matched_layer.kappa.z[(K0.index)] * matched_layer.gamma[k] / (2. * params.dt) * C.b.z[K0.index] *
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
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
              if (is_conductive)
                J_c.yz[(K1.index)][j][i] -=
                        rho_cond.z[(K1.index)] * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                              IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
              if (params.is_disp_ml)
                J_s.yz[(K1.index)][j][i] +=
                        matched_layer.kappa.z[(K1.index)] * matched_layer.gamma[k] / (2. * params.dt) * C.b.z[K1.index] *
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
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
              if (is_conductive)
                J_c.xz[(K1.index)][j][i] +=
                        rho_cond.z[(K1.index)] * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                              IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
              if (params.is_disp_ml)
                J_s.xz[(K1.index)][j][i] -=
                        matched_layer.kappa.z[(K1.index)] * matched_layer.gamma[k] / (2. * params.dt) * C.b.z[K1.index] *
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
                               exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
                          exp(-1.0 * DCPI *
                              pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
          //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
          if (is_conductive)
            J_c.yz[K0.index][j][i] +=
                    rho_cond.z[K0.index] * C.b.z[K0.index] *
                    real((Ksource.real[0][i - (I0.index)][2] +
                          IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][2]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
                    exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
          //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
          if (params.is_disp_ml) {
            J_s.yz[K0.index][j][i] -=
                    matched_layer.kappa.z[K0.index] * matched_layer.gamma[K0.index] / (2. * params.dt) *
                    C.b.z[K0.index] *
                    real((Ksource.real[0][i - (I0.index)][2] +
                          IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][2]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
                    exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
            //J_s.yz[(int)K0[0]][j][i] -= matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
          }
        }
      } else
        for (j = 0; j < J_tot; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            /*
        if(i==41 & j==41)
        fprintf(stderr,"C.b.z = %.10e, Re(K) = %.10e, Im(K) = %.10e, iteration_controller.time_H= %.10e, params.to_l=%.10e, params.delta.dz/LIGHT_V/2=%.10e, hwhm = %.10e, dE=%.10e\n",C.b.z[(int)K0[0]],Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2],Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2],iteration_controller.time_H,params.to_l,params.delta.dz/LIGHT_V/2,params.hwhm,C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l + params.delta.dz/LIGHT_V/2.)/(params.hwhm),2)));
      */
            E_s.yz[K0.index][j][i] =
                    E_s.yz[K0.index][j][i] -
                    C.b.z[K0.index] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                                  IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                                 (-1.0 * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
                            exp(-1.0 * DCPI *
                                pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
            //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
            if (is_conductive)
              J_c.yz[K0.index][j][i] +=
                      rho_cond.z[K0.index] * C.b.z[K0.index] *
                      real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                            IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                           (-1.0 * IMAGINARY_UNIT) *
                           exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
                      exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
            //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
            if (params.is_disp_ml) {
              J_s.yz[K0.index][j][i] -=
                      matched_layer.kappa.z[K0.index] * matched_layer.gamma[K0.index] / (2. * params.dt) *
                      C.b.z[K0.index] *
                      real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                            IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                           (-1.0 * IMAGINARY_UNIT) *
                           exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
                      exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
              //J_s.yz[(int)K0[0]][j][i] -= matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
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
                               exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2 * DCPI))) *
                          exp(-1.0 * DCPI *
                              pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
          //E_s.xz[(int)K0[0]][j][i] = E_s.xz[(int)K0[0]][j][i] + C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2 ));
          if (is_conductive)
            J_c.xz[K0.index][j][i] -=
                    rho_cond.z[K0.index] * C.b.z[K0.index] *
                    real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                          IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2 * DCPI))) *
                    exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
          //J_c.xz[(int)K0[0]][j][i] -= rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2 ));
          if (params.is_disp_ml)
            J_s.xz[K0.index][j][i] +=
                    matched_layer.kappa.z[K0.index] * matched_layer.gamma[K0.index] / (2. * params.dt) *
                    C.b.z[K0.index] *
                    real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                          IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2 * DCPI))) *
                    exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
          //J_s.xz[(int)K0[0]][j][i] += matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2 ));
        }
      //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
      H.ft = real((-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_H - params.to_l), 2. * DCPI))) *
             exp(-1.0 * DCPI * pow((iteration_controller.time_H - params.to_l + params.delta.dz / LIGHT_V / 2.) / (params.hwhm), 2));
      //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(iteration_controller.time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((iteration_controller.time_H - params.to_l)/(params.hwhm),2));
    }
    //fprintf(stderr,"Pos 10:\n");

    //end of source terms
    if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }

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
    fprintf(stdout,"iteration_controller.tind: %d\n",iteration_controller.tind);
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
                  /*if(iteration_controller.tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                  H_s.yz[k][j][i] = D.a.z[k_loc] * H_s.yz[k][j][i] +
                                    D.b.z[k_loc] * (E_s.xy[k][j][i] + E_s.xz[k][j][i] -
                                                    E_s.xy[k + 1][j][i] - E_s.xz[k + 1][j][i]);
                } else {
                  /*if(iteration_controller.tind==0)
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
                  /*if(iteration_controller.tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                  //H_s.yz[k][j][i] = D.a.z[k_loc]*H_s.yz[k][j][i]+D.b.z[k_loc]*(E_s.xy[k][j][i] + E_s.xz[k][j][i] - E_s.xy[k+1][j][i] - E_s.xz[k+1][j][i]);
                } else {
                  ca_vec[n][k] = Dmaterial.a.z[materials[k][j][i] - 1];
                  cb_vec[n][k] = Dmaterial.b.z[materials[k][j][i] - 1];
                  /*if(iteration_controller.tind==0)
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
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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


            first_derivative(eh_vec[n], eh_vec[n], dk_h_x,
                             N_h_x, H_s.zx.plan_f[n], H_s.zx.plan_b[n]);

            for (i = 0; i < I_tot; i++) {
              H_s.zx[k][j][i] = ca_vec[n][i] * H_s.zx[k][j][i] -
                                cb_vec[n][i] *
                                        eh_vec[n][i][0] / ((double) N_h_x);
            }
          }
//PSTD, H_s.zx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
      }//(params.dimension==THREE || params.dimension==TE)
    }  //end parallel
    if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }

    //fprintf(stderr,"Pos 11b:\n");
    //update terms for self consistency across scattered/total interface - E updates
    if (params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase = exp(-IMAGINARY_UNIT * fmod(params.omega_an * iteration_controller.time_E, 2. * DCPI));
      double commonAmplitude = linearRamp(iteration_controller.time_E, 1. / (params.omega_an / (2 * DCPI)), ramp_width);
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
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][0]));
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              H_s.yx[k][j][(I0.index) - 1] =
                      H_s.yx[k][j][(I0.index) - 1] -
                      D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][1] +
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][1]));
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
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][4]));
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              H_s.yx[k][j][(I1.index)] =
                      H_s.yx[k][j][(I1.index)] +
                      D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][5] +
                                    IMAGINARY_UNIT * Isource.imag[k - (K0.index)][j - (J0.index)][5]));
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
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][0]));

            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              H_s.xy[k][(J0.index) - 1][i] =
                      H_s.xy[k][(J0.index) - 1][i] +
                      D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][1] +
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][1]));
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
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][4]));
            if (k < (K1.index) || params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              H_s.xy[k][(J1.index)][i] =
                      H_s.xy[k][(J1.index)][i] -
                      D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][5] +
                                    IMAGINARY_UNIT * Jsource.imag[k - (K0.index)][i - (I0.index)][5]));
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
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][0]));
            if (j < (J1.index))
              H_s.xz[(K0.index) - 1][j][i] =
                      H_s.xz[(K0.index) - 1][j][i] -
                      D.b.z[(K0.index) - 1] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][1] +
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][1]));
          }
          if (K1.apply) {//Perform across K1
            if (i < (I1.index))
              H_s.yz[(K1.index)][j][i] =
                      H_s.yz[(K1.index)][j][i] -
                      D.b.z[(K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][4] +
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][4]));
            if (j < (J1.index))
              H_s.xz[(K1.index)][j][i] =
                      H_s.xz[(K1.index)][j][i] +
                      D.b.z[(K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][5] +
                                    IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][5]));
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
                               exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_E - params.to_l), 2 * DCPI))) *
                          exp(-1. * DCPI * pow((iteration_controller.time_E - params.to_l) / (params.hwhm), 2.));
          //broadband source term
          if (params.eyi_present)
            H_s.xz[(K0.index) - 1][j][i] =
                    H_s.xz[(K0.index) - 1][j][i] - D.b.z[(K0.index) - 1] * Ei.y[iteration_controller.tind][j][i];
        }
        //fprintf(stderr,"Pos 11e\n");
        for (i = 0; i < I_tot; i++) {
          H_s.yz[(K0.index) - 1][j][i] =
                  H_s.yz[(K0.index) - 1][j][i] +
                  D.b.z[(K0.index) - 1] *
                          real((Ksource.real[0][i - (I0.index)][0] +
                                IMAGINARY_UNIT * Ksource.imag[0][i - (I0.index)][0]) *
                               (-1. * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_E - params.to_l), 2 * DCPI))) *
                          exp(-1. * DCPI * pow((iteration_controller.time_E - params.to_l) / (params.hwhm), 2.));
          //broadband source term
          if (params.exi_present)
            H_s.yz[(K0.index) - 1][j][i] =
                    H_s.yz[(K0.index) - 1][j][i] + D.b.z[(K0.index) - 1] * Ei.x[iteration_controller.tind][j][i];
          //if(i==511)
          //  fprintf(stdout,"%e\n",D.b.z[((int)K0[0])-1]*exi[iteration_controller.tind][j][i]);
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
                                  IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][1]) *
                                 (-1. * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_E - params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((iteration_controller.time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.eyi_present)
              H_s.xz[(K0.index) - 1][j][i] =
                      H_s.xz[(K0.index) - 1][j][i] - D.b.z[(K0.index) - 1] * Ei.y[iteration_controller.tind][j][i];
          }
        //fprintf(stderr,"Pos 11h\n");
        for (j = 0; j < (J_tot + 1); j++)
          for (i = 0; i < I_tot; i++) {
            H_s.yz[(K0.index) - 1][j][i] =
                    H_s.yz[(K0.index) - 1][j][i] +
                    D.b.z[(K0.index) - 1] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][0] +
                                  IMAGINARY_UNIT * Ksource.imag[j - (J0.index)][i - (I0.index)][0]) *
                                 (-1. * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_E - params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((iteration_controller.time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.exi_present)
              H_s.yz[(K0.index) - 1][j][i] =
                      H_s.yz[(K0.index) - 1][j][i] + D.b.z[(K0.index) - 1] * Ei.x[iteration_controller.tind][j][i];
          }
        //fprintf(stderr,"Pos 11i\n");
      }
      E.ft = real((-1. * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(params.omega_an * (iteration_controller.time_E - params.to_l), 2 * DCPI))) *
             exp(-1. * DCPI * pow((iteration_controller.time_E - params.to_l) / (params.hwhm), 2.));
      //fprintf(stderr,"Pos 11j\n");
    }
    if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL);; }

    if (params.exphasorssurface || params.exphasorsvolume || params.exdetintegral || vertex_phasors.there_are_vertices_to_extract_at()) {
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

      if( (iteration_controller.tind % iteration_controller.Nsteps) == 0 ){
        E.angular_norm = 0.0;
        H.angular_norm = 0.0;

        for(int ifx=0; ifx<f_ex_vec.size(); ifx++){
          E_norm[ifx] = 0.;
          H_norm[ifx] = 0.;
        }
      }

	/*In the calls below, the following two lines of code are equivalent up to numerical precision:

	  E.add_to_angular_norm(fte, iteration_controller.tind, iteration_controller.Nsteps, params);
	  E.add_to_angular_norm(fte, iteration_controller.tind % iteration_controller.Nsteps, iteration_controller.Nsteps, params);

	  To understand why, first consult the lines of code above:

	  Nsteps_tmp = ceil(2.*DCPI/omega_an[0]/dt[0]*3);
	  dt[0] = 2.*DCPI/omega_an[0]*3/Nsteps_tmp;
	  iteration_controller.Nsteps = (int)lround(Nsteps_tmp);

	  Where dt and iteration_controller.Nsteps are set. The reason for the factor of 3 is that we will perform complex
	  amplitude extraction over 3 fundamental periods of the monochromatic source. We can then make
	  the following statement:

	  T/dt*3=1/(f*dt)*3=iteration_controller.Nsteps

	  where T and f (omega=2*pi*f) are the period and frequency of the monochromatic source, respectively.

	  Then consider the argument of the exponentional function on phasor_norm, called by add_to_angular_norm, where iteration_controller.tind=n is used:

	  i*omega*((double) (n+1))*dt (where fmod(.,2*DCPI) is ignored since this will not affect the result)

	  The argument of this function simplifies to:

	  i*omega*(iteration_controller.tind+1)*dt=i*2*pi*f*(iteration_controller.tind+1)*dt=i*2*pi*(iteration_controller.tind+1)*3/iteration_controller.Nsteps (using f*dt=3/iteration_controller.Nsteps)

	  Then, without loss of generallity, let iteration_controller.tind = p*iteration_controller.Nsteps + q, substituting into the above

	  i*2*pi*(iteration_controller.tind+1)*3/iteration_controller.Nsteps = i*2*pi*(p*iteration_controller.Nsteps + q)*3/iteration_controller.Nsteps = i*2*pi*3*p + i*2*pi*q*3/iteration_controller.Nsteps

	  In which case exp(i*2*pi*3*p + i*2*pi*q*3/iteration_controller.Nsteps) = exp(i*2*pi*q*3/iteration_controller.Nsteps)

	  If instead we use iteration_controller.tind % iteration_controller.Nsteps=n, we see that n=q, leading to the same exponential function as
	  above. So the two cases are equivalent.
	 */

        E.add_to_angular_norm(iteration_controller.tind, iteration_controller.Nsteps, params);
        H.add_to_angular_norm(iteration_controller.tind, iteration_controller.Nsteps, params);

        for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
          extractPhasorENorm(&E_norm[ifx], E.ft, iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI, params.dt, iteration_controller.Nsteps);
          extractPhasorHNorm(&H_norm[ifx], H.ft, iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI, params.dt, iteration_controller.Nsteps);
        }
      } else {
        if ((iteration_controller.tind - params.start_tind) % params.Np == 0) {

          E.add_to_angular_norm(iteration_controller.tind, params.Npe, params);
          H.add_to_angular_norm(iteration_controller.tind, params.Npe, params);

          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            extractPhasorENorm(&E_norm[ifx], E.ft, iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI, params.dt, params.Npe);
            extractPhasorHNorm(&H_norm[ifx], H.ft, iteration_controller.tind, f_ex_vec[ifx] * 2 * DCPI, params.dt, params.Npe);
          }
        }
      }
    }
    if (TIME_EXEC) { iteration_controller.click_timer(IterationTimers::INTERNAL); }

    iteration_controller.log_update(E_s, H_s, iteration_controller.tind);
    //fprintf(stderr,"Post-iter 3\n");
    if ((params.source_mode == SourceMode::steadystate) && (iteration_controller.tind == (params.Nt - 1)) && (params.run_mode == RunMode::complete) &&
        params.exphasorsvolume) {
      fprintf(stdout, "Iteration limit reached, setting output fields to last complete DFT\n");
      E.set_values_from(E_copy);
    }
    //fprintf(stderr,"Post-iter 4\n");
    fflush(stdout);
    //fprintf(stderr,"Post-iter 5\n");
    //fprintf(stderr,"%s %d %d\n",tdfdirstr, strcmp(tdfdirstr,""),are_equal(tdfdirstr,""));
    if (params.has_tdfdir && (iteration_controller.tind % params.Np) == 0) {
      fprintf(stderr,"Saving field\n");
      ex_td_field_exporter.export_field(E_s, skip_tdf, iteration_controller.tind);
    }
    //fprintf(stderr,"Post-iter 6\n");
    /*write out fdtdgrid to a file*/
    /*
     MATFile *toutfile;
     char toutputfilename[100];
     if(iteration_controller.tind % params.Np == 0){
     //if(iteration_controller.tind <= 1000){
       sprintf(toutputfilename,"tdata/fdtdgrid_%04d.mat",tind);
       toutfile = matOpen(toutputfilename, "w");
       matPutVariable(toutfile, "fdtdgrid", (mxArray *)in_matrices[0]);
       matClose(toutfile);
       }
    */
    /*write out fdtdgrid to a file*/

  }//end of main iteration loop
  if (TIME_MAIN_LOOP) {
    //fprintf(stderr,"Post-iter 7\n");
    iteration_controller.end_timer(IterationTimers::MAIN);
    //fprintf(stderr,"Post-iter 8\n");
    spdlog::info("Time (seconds) ellapsed in main loop: {0:e}", iteration_controller.time_ellapsed_by(IterationTimers::MAIN));
    //fprintf(stderr,"Post-iter 9\n");
  }
  //save state of fdtdgrid

  //fprintf(stderr,"Pos 12\n");
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    E.normalise_volume();
    H.normalise_volume();
  }

  //fprintf(stderr,"Pos 13\n");
  if (params.run_mode == RunMode::complete && params.exphasorssurface) {
    spdlog::info("Surface phasors");
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      surface_phasors.normalise_surface(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
    }
  }
  if (params.run_mode == RunMode::complete && vertex_phasors.there_are_vertices_to_extract_at()) {
    spdlog::info("Vertex phasors");
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      vertex_phasors.normalise_vertices(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
    }
  }

  //fprintf(stderr,"Pos 14\n");
  if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete && params.exdetintegral) {
    for (int im = 0; im < D_tilde.num_det_modes(); im++)
      for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
        Idx[ifx][im] = Idx[ifx][im] / E_norm[ifx];
        Idy[ifx][im] = Idy[ifx][im] / E_norm[ifx];

        Idx_re[ifx][im] = real(Idx[ifx][im]);
        Idx_im[ifx][im] = imag(Idx[ifx][im]);

        Idy_re[ifx][im] = real(Idy[ifx][im]);
        Idy_im[ifx][im] = imag(Idy[ifx][im]);
      }
  }

  //now find the maximum absolute value of residual field in the grid
  // after resetting the maxfield value calculated in the main loop
  maxfield = max(E_s.largest_field_value(), H_s.largest_field_value());
  //fprintf(stderr,"Pos 15\n");
  //noe set the output
  ndims = 2;
  dims[0] = 1;
  dims[1] = 1;
  plhs[25] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  *mxGetPr((mxArray *) plhs[25]) = maxfield;

  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    output_grid_labels.initialise_from(input_grid_labels, E.il, E.iu, E.jl, E.ju, E.kl, E.ku);
  }

  auto interp_output_grid_labels = GridLabels();

  //fprintf(stderr,"Pos 15_m1\n");
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    //now interpolate over the extracted phasors
    if (params.dimension == THREE) {
      E.interpolate_over_range(&plhs[13], &plhs[14], &plhs[15], 2, E.I_tot - 2, 2, E.J_tot - 2, 2,
                               E.K_tot - 2, Dimension::THREE);
      H.interpolate_over_range(&plhs[16], &plhs[17], &plhs[18], 2, H.I_tot - 2, 2, H.J_tot - 2, 2,
                               H.K_tot - 2, Dimension::THREE);
    } else {
      // either TE or TM, but interpolate_over_range will handle that for us. Only difference is the k_upper/lower values we pass...
      E.interpolate_over_range(&plhs[13], &plhs[14], &plhs[15], 2, E.I_tot - 2, 2, E.J_tot - 2, 0,
                               0, params.dimension);
      H.interpolate_over_range(&plhs[16], &plhs[17], &plhs[18], 2, H.I_tot - 2, 2, H.J_tot - 2, 0,
                               0, params.dimension);
    }

    //fprintf(stderr,"Pos 15a\n");
    //now set up the grid labels for the interpolated fields
    label_dims[0] = 1;
    label_dims[1] = E.I_tot - 3;
    plhs[19] = mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//x
    //fprintf(stderr,"Pos 15b\n");
    label_dims[0] = 1;
    label_dims[1] = E.J_tot - 3;
    if (label_dims[1] < 1) label_dims[1] = 1;
    //fprintf(stderr,"creating plhs[20]: %d,%d\n",label_dims[0],label_dims[1]);
    plhs[20] = mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//y
    //fprintf(stderr,"Pos 15c\n");
    label_dims[0] = 1;
    if (params.dimension == THREE) label_dims[1] = E.K_tot - 3;
    else
      label_dims[1] = 1;
    //fprintf(stderr,"Pos 15d\n");
    plhs[21] = mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//z
    //fprintf(stderr,"Pos 15e\n");

    interp_output_grid_labels.x = mxGetPr((mxArray *) plhs[19]);
    interp_output_grid_labels.y = mxGetPr((mxArray *) plhs[20]);
    interp_output_grid_labels.z = mxGetPr((mxArray *) plhs[21]);

    if (params.dimension == THREE) {
      interp_output_grid_labels.initialise_from(output_grid_labels, 2, E.I_tot - 2, 2, E.J_tot - 2,
                                                2, E.K_tot - 2);
    } else {
      interp_output_grid_labels.initialise_from(output_grid_labels, 2, E.I_tot - 2, 2, E.J_tot - 2,
                                                0, 0);
    }
    //fprintf(stderr,"Pos 15f\n");
  } else {
    mwSize *emptydims;
    emptydims = (mwSize *) malloc(2 * sizeof(mwSize));
    int emptyloop;
    emptydims[0] = 0;
    emptydims[1] = 0;
    for (emptyloop = 13; emptyloop <= 18; emptyloop++)
      plhs[emptyloop] =
              mxCreateNumericArray(2, (const mwSize *) emptydims, mxDOUBLE_CLASS, mxCOMPLEX);
    for (emptyloop = 19; emptyloop <= 21; emptyloop++)
      plhs[emptyloop] =
              mxCreateNumericArray(2, (const mwSize *) emptydims, mxDOUBLE_CLASS, mxCOMPLEX);
    free(emptydims);
  }


  //fprintf(stderr,"Pos 16\n");
  /*Now export 3 matrices, a vertex list, a matrix of complex amplitudes at
    these vertices and a list of facets*/
  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    //first regenerate the mesh since we threw away the facet list before iterating
    mxArray *dummy_vertex_list;
    if (J_tot == 0)
      conciseCreateBoundary(cuboid[0], cuboid[1], cuboid[4], cuboid[5], &dummy_vertex_list,
                            &mx_surface_facets);
    else
      conciseTriangulateCuboidSkip(cuboid[0], cuboid[1], cuboid[2], cuboid[3], cuboid[4], cuboid[5],
                                   params.spacing_stride, &dummy_vertex_list,
                                   &mx_surface_facets);
    mxDestroyArray(dummy_vertex_list);

    //now create and populate the vertex list
    surface_phasors.create_vertex_list(input_grid_labels);

    //assign outputs
    plhs[22] = surface_phasors.get_vertex_list();
    plhs[23] = surface_phasors.get_mx_surface_amplitudes();
    plhs[24] = mx_surface_facets;

  } else {//still set outputs
    ndims = 2;
    dims[0] = 0;
    dims[1] = 0;
    plhs[22] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
    plhs[23] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
    plhs[24] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  }


  /*End of FDTD iteration*/

  //fprintf(stderr,"Pos 17\n");
  /* Free the additional data structures used to cast the matlab arrays*/
  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    mxDestroyArray(mx_surface_vertices);
    // ~SurfacePhasors cleans up remaining surface phasor arrays, except the data which we return in plhs[23]
  }
  // ~VertexPhasors cleans up the vertex phasors arrays, except the data which we return in plhs[28]

  // some smart AF cleanup is going to be needed here...
  if (params.exdetintegral) {
    free_cast_matlab_2D_array(Idx_re);
    free_cast_matlab_2D_array(Idx_im);
    free_cast_matlab_2D_array(Idy_re);
    free_cast_matlab_2D_array(Idy_im);
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      free(Idx[ifx]);
      free(Idy[ifx]);
    }
    free(Idx);
    free(Idy);
    /*
      for(int j=0;j<(J_tot-params.pml.Dyl-params.pml.Dyu);j++){
      free(Ex_t.cm[j]);free(Ey_t.cm[j]);
      }
      //fprintf(stderr,"Position 9\n");
      free(Ex_t_cm);free(Ey_t_cm);
      //fprintf(stderr,"Position 10\n");
    */
  }

  //fprintf(stderr,"Pos 20\n");

  free_cast_matlab_2D_array(iwave_lEx_Rbs);
  free_cast_matlab_2D_array(iwave_lEx_Ibs);
  free_cast_matlab_2D_array(iwave_lEy_Rbs);
  free_cast_matlab_2D_array(iwave_lEy_Ibs);

  free_cast_matlab_2D_array(iwave_lHx_Rbs);
  free_cast_matlab_2D_array(iwave_lHx_Ibs);
  free_cast_matlab_2D_array(iwave_lHy_Rbs);
  free_cast_matlab_2D_array(iwave_lHy_Ibs);

  if (params.dimension == THREE) free_cast_matlab_3D_array(materials, E_s.K_tot + 1);
  else
    free_cast_matlab_3D_array(materials, 0);
    /*Free the additional memory which was allocated to store integers which were passed as doubles*/


  //must destroy surface_phasors.mx_surface_amplitudes
}

/*Sets the contents of the 3 dimensional double array to zero
  inArray - pointer to the array
  i_lim - number of elements along the i dimension of the array
  j_lim - number of elements along the j dimension of the array
  k_lim - number of elements along the k dimension of the array

  The array is assumed to be indexed according to inArray[k][j][i]

*/

void initialiseDouble3DArray(double ***inArray, int i_lim, int j_lim, int k_lim) {
  for (int k_var = 0; k_var < k_lim; k_var++)
    for (int j_var = 0; j_var < j_lim; j_var++)
      for (int i_var = 0; i_var < i_lim; i_var++) inArray[k_var][j_var][i_var] = 0.0;
}

/*Sets the contents of the 2 dimensional double array to zero
  inArray - pointer to the array
  i_lim - number of elements along the i dimension of the array
  j_lim - number of elements along the j dimension of the array

  The array is assumed to be indexed according to inArray[j][i]

*/

void initialiseDouble2DArray(double **inArray, int i_lim, int j_lim) {
  for (int j_var = 0; j_var < j_lim; j_var++)
    for (int i_var = 0; i_var < i_lim; i_var++) inArray[j_var][i_var] = 0.0;
}

void extractPhasorENorm(complex<double> *Enorm, double ft, int n, double omega, double dt, int Nt) {
  *Enorm += ft * exp(fmod(omega * ((double) (n + 1)) * dt, 2 * DCPI) * IMAGINARY_UNIT) * 1. / ((double) Nt);
}

void extractPhasorHNorm(complex<double> *Hnorm, double ft, int n, double omega, double dt, int Nt) {
  *Hnorm += ft * exp(fmod(omega * ((double) n + 0.5) * dt, 2 * DCPI) * IMAGINARY_UNIT) * 1. / ((double) Nt);
}

void extractPhasorsPlane(double **iwave_lEx_Rbs, double **iwave_lEx_Ibs, double **iwave_lEy_Rbs,
                         double **iwave_lEy_Ibs, double **iwave_lHx_Rbs, double **iwave_lHx_Ibs,
                         double **iwave_lHy_Rbs, double **iwave_lHy_Ibs,
                         ElectricSplitField &E, MagneticSplitField &H, int I_tot, int J_tot, int K1, int n,
                         double omega, double dt, int Nt) {

  complex<double> phaseTerm = 0., subResult = 0.;


  phaseTerm = fmod(omega * ((double) n) * dt, 2 * DCPI);
  int i, j;

  for (j = 0; j < J_tot; j++)
    for (i = 0; i < (I_tot + 1); i++) {


      //Eyz
      subResult = (E.yz[K1][j][i] + E.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);

      iwave_lEy_Rbs[j][i] = iwave_lEy_Rbs[j][i] + real(subResult);
      iwave_lEy_Ibs[j][i] = iwave_lEy_Ibs[j][i] + imag(subResult);

      //Hxz
      subResult = (H.xz[K1 - 1][j][i] + H.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);

      iwave_lHx_Rbs[j][i] = iwave_lHx_Rbs[j][i] + real(subResult);
      iwave_lHx_Ibs[j][i] = iwave_lHx_Ibs[j][i] + imag(subResult);
    }

  for (j = 0; j < (J_tot + 1); j++)
    for (i = 0; i < I_tot; i++) {


      //Exz
      subResult = (E.xz[K1][j][i] + E.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);

      iwave_lEx_Rbs[j][i] = iwave_lEx_Rbs[j][i] + real(subResult);
      iwave_lEx_Ibs[j][i] = iwave_lEx_Ibs[j][i] + imag(subResult);

      //Hyz
      subResult = (H.yz[K1 - 1][j][i] + H.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);

      iwave_lHy_Rbs[j][i] = iwave_lHy_Rbs[j][i] + real(subResult);
      iwave_lHy_Ibs[j][i] = iwave_lHy_Ibs[j][i] + imag(subResult);
    }
}

/*Implements a linear ramp which has the properties:

  ramp(t) = 1 if t > rampwidth*period
  = t/(rampwidth*period) otherwise

  t - the current time at which to evaluate the ramp
  period - the period of the monochormatic sinusoidal excitation
  rampwidth - the fraction or number of periods with of the ramp

*/
double linearRamp(double t, double period, double rampwidth) {

  if (t > period * rampwidth) return 1.;
  else
    return t / (period * rampwidth);
}

/* These functions are used by the dispersive component of the code*/

/*Work out if there are any non-zero alpha values*/
bool is_dispersive(unsigned char ***materials, double *gamma, double dt, int I_tot, int J_tot,
                  int K_tot) {
  int max_mat = 0;

  //first find the number of entries in alpha
  for (int k = 0; k < (K_tot + 1); k++)
    for (int j = 0; j < (J_tot + 1); j++)
      for (int i = 0; i < (I_tot + 1); i++) {
        if (materials[k][j][i] > max_mat) max_mat = materials[k][j][i];
      }
  //now see if there are any non zero alphas
  for (int i = 0; i < max_mat; i++) {
    if (fabs(gamma[i] / dt) > 1e-15) { return 1; }
  }
  return 0;
}
