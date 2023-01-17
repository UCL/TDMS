#include "iterator.h"

#include <complex>
#include <algorithm>
#include <ctime>

#include <omp.h>
#include <spdlog/spdlog.h>

#include "mesh_base.h"
#include "numerical_derivative.h"
#include "array_init.h"
#include "fieldsample.h"
#include "globals.h"
#include "interface.h"
#include "matlabio.h"
#include "shapes.h"
#include "source.h"
#include "surface_phasors.h"
#include "vertex_phasors.h"
#include "objects_from_infile.h"
#include "pstd_variables.h"
#include "output_matrices.h"
#include "loop_timers.h"
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
OutputMatrices execute_simulation(InputMatrices in_matrices, SolverMethod solver_method,
                                  PreferredInterpolationMethods preferred_interpolation_methods) {
  // declare the to-be output
  OutputMatrices outputs;

  // PSTD storage (to be a member of the simulation manager later)
  PSTDVariables PSTD;

  // timers for the main loop (to be a member of the simulation manager later)
  LoopTimers timers;

  // this will later be a member of simulation_elements
  ObjectsFromInfile inputs(in_matrices, solver_method, preferred_interpolation_methods);

  spdlog::info("Using {} OMP threads", omp_get_max_threads());

  auto J_s = CurrentDensitySplitField();

  auto E_copy = ElectricField();  // Used to check convergence with E - E_copy
  E_copy.set_preferred_interpolation_methods(preferred_interpolation_methods); // We never actually interpolate this field, but adding this just in case we later add functionality that depends on it

  double rho;
  double alpha_l, beta_l, gamma_l;
  double kappa_l, sigma_l;
  double t0;

  double Ca, Cb, Cc;     //used by interpolation scheme
  //the C and D vars for free space and pml
  double Enp1, Jnp1;

  double maxfield = 0.;

  //refractive index of the first layer of the multilayer, or of the bulk of homogeneous
  double refind;
  double phaseTermE;
  complex<double> cphaseTermE;
  double lambda_an_t;

  int i, j, k, n, k_loc, K;
  int Nsteps = 0, dft_counter = 0;

  mxArray *E_copy_MATLAB_data[3];
  mxArray *mx_surface_vertices, *mx_surface_facets;
  complex<double> Idxt, Idyt, kprop;

  // get the number of Yee cells in each axial direction
  IJKDims IJK_tot = inputs.IJK_tot;
  outputs.set_n_Yee_cells(IJK_tot);
  int I_tot = IJK_tot.I_tot(), J_tot = IJK_tot.J_tot(), K_tot = IJK_tot.K_tot();

  /*Deduce the refractive index of the first layer of the multilayer, or of the bulk of homogeneous*/
  refind = sqrt(1. / (inputs.freespace_Cbx[0] / inputs.params.dt * inputs.params.delta.dx) / EPSILON0);
  spdlog::info("refind = {0:e}", refind);

  /*Setup temporary storage for detector sensitivity evaluation*/
  auto Ex_t = DetectorSensitivityArrays();
  auto Ey_t = DetectorSensitivityArrays();

  if (inputs.params.exdetintegral) {
    int n0 = I_tot - inputs.params.pml.Dxl - inputs.params.pml.Dxu;
    int n1 = J_tot - inputs.params.pml.Dyl - inputs.params.pml.Dyu;
    Ex_t.initialise(n1, n0);
    Ey_t.initialise(n1, n0);
  }

  // setup PSTD variables, and any dependencies there might be
  PSTD.set_using_dimensions(IJK_tot);
  EHVec eh_vec;
  if (solver_method == SolverMethod::PseudoSpectral) {
    int max_IJK = IJK_tot.max_IJK(), n_threads = omp_get_max_threads();
    eh_vec.allocate(n_threads, max_IJK + 1);

    inputs.E_s.initialise_fftw_plan(n_threads, eh_vec);
    inputs.H_s.initialise_fftw_plan(n_threads, eh_vec);
  }

  //initialise E_norm and H_norm
  auto E_norm = (complex<double> *) malloc(inputs.f_ex_vec.size() * sizeof(complex<double>));
  auto H_norm = (complex<double> *) malloc(inputs.f_ex_vec.size() * sizeof(complex<double>));
  for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
    E_norm[ifx] = 0.;
    H_norm[ifx] = 0.;
  }

  //fprintf(stderr,"Pre 02\n");

  //fprintf(stderr,"Qos 00 (%d) (%d,%d,%d,%d):\n",J_tot,cuboid[0], cuboid[1],cuboid[4], cuboid[5]);
  /*set up surface mesh if required*/

  if (inputs.params.exphasorssurface && inputs.params.run_mode == RunMode::complete) {
    if (J_tot == 0)
      conciseCreateBoundary(inputs.cuboid[0], inputs.cuboid[1], inputs.cuboid[4], inputs.cuboid[5], &mx_surface_vertices,
                            &mx_surface_facets);
    else
      conciseTriangulateCuboidSkip(inputs.cuboid[0], inputs.cuboid[1], inputs.cuboid[2], inputs.cuboid[3], inputs.cuboid[4], inputs.cuboid[5],
                                   inputs.params.spacing_stride, &mx_surface_vertices,
                                   &mx_surface_facets);
    //fprintf(stderr,"Qos 00a:\n");
    //we don't need the facets so destroy the matrix now to save memory
    mxDestroyArray(mx_surface_facets);

    outputs.surface_phasors.set_from_matlab_array(mx_surface_vertices, inputs.f_ex_vec.size());
  }

  //fprintf(stderr,"Pre 03\n");
  /*Now set up the phasor array, we will have 3 complex output arrays for Ex, Ey and Ez.
    Phasors are extracted over the range Dxl + 3 - 1 to I_tot - Dxu - 1 to avoid pml cells
    see page III.80 for explanation of the following. This has been extended so that interpolation
    is done at the end of the FDTD run and also to handle the case of when there is no PML in place
    more appropriatley*/
  outputs.setup_EH_and_gridlabels(inputs.params, inputs.input_grid_labels, preferred_interpolation_methods);

  // We need to test for convergence under the following conditions. As such, we need to initialise the array that will ultimately be copies of the phasors at the previous iteration, to test convergence against
  if (inputs.params.run_mode == RunMode::complete && inputs.params.exphasorsvolume &&
      inputs.params.source_mode == SourceMode::steadystate) {
    // fetch the E-field dimensions from the output object
    IJKDims IJK = outputs.get_E_dimensions();
    int dummy_dims[3] = {IJK.I_tot(), IJK.J_tot(), IJK.K_tot()};
    // allocate memory space for this array
    E_copy_MATLAB_data[0] =
            mxCreateNumericArray(3, (const mwSize *) dummy_dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
    E_copy_MATLAB_data[1] =
            mxCreateNumericArray(3, (const mwSize *) dummy_dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
    E_copy_MATLAB_data[2] =
            mxCreateNumericArray(3, (const mwSize *) dummy_dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez

    E_copy.real.x = cast_matlab_3D_array(mxGetPr(E_copy_MATLAB_data[0]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);
    E_copy.imag.x = cast_matlab_3D_array(mxGetPi(E_copy_MATLAB_data[0]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);

    E_copy.real.y = cast_matlab_3D_array(mxGetPr(E_copy_MATLAB_data[1]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);
    E_copy.imag.y = cast_matlab_3D_array(mxGetPi(E_copy_MATLAB_data[1]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);

    E_copy.real.z = cast_matlab_3D_array(mxGetPr(E_copy_MATLAB_data[2]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);
    E_copy.imag.z = cast_matlab_3D_array(mxGetPi(E_copy_MATLAB_data[2]), dummy_dims[0],
                                         dummy_dims[1], dummy_dims[2]);

    E_copy.I_tot = IJK.I_tot();
    E_copy.J_tot = IJK.J_tot();
    E_copy.K_tot = IJK.K_tot();

    E_copy.zero();
  }

  // Setup the ID output
  bool need_Id_memory = (inputs.params.exdetintegral && inputs.params.run_mode == RunMode::complete);
  outputs.setup_Id(!need_Id_memory, inputs.f_ex_vec.size(), inputs.D_tilde.num_det_modes());

  /*This is just for efficiency */
  K = K_tot - inputs.params.pml.Dxl - inputs.params.pml.Dxu;

  /*Now set up the phasor arrays for storing the fdtd version of the input fields,
    these will be used in a boot strapping procedure. Calculated over a complete
    xy-plane. */

  // these are needed later... but don't seem to EVER be used? They were previously plhs[6->9], but these outputs were never written. Also, they are assigned to, but never written out nor referrenced by any of the other variables in the main loop. I am confused... Also note that because we're using the Matrix class, we order indices [i][j][k] rather than [k][j][i] like in the rest of the codebase :(
  Matrix<complex<double>> iwave_lEx_bs, iwave_lEy_bs, iwave_lHx_bs, iwave_lHy_bs;

  //x electric field source phasor - boot strapping
  iwave_lEx_bs.allocate(I_tot, J_tot + 1);
  //y electric field source phasor - boot strapping
  iwave_lEy_bs.allocate(I_tot + 1, J_tot);
  //x magnetic field source phasor - boot strapping
  iwave_lHx_bs.allocate(I_tot + 1, J_tot);
  //y magnetic field source phasor - boot strapping
  iwave_lHy_bs.allocate(I_tot, J_tot + 1);

  /*start dispersive*/

  //work out if we have any disperive materials
  bool is_disp = is_dispersive(inputs.materials, inputs.gamma, inputs.params.dt, I_tot, J_tot, K_tot);
  //work out if we have conductive background: background is conductive if at least one entry exceeds 1e-15
  bool is_conductive = !(inputs.rho_cond.all_elements_less_than(1e-15, I_tot + 1, J_tot + 1, K_tot + 1));
  // work out if we have a dispersive background
  if (inputs.params.is_disp_ml) inputs.params.is_disp_ml = inputs.matched_layer.is_dispersive(K_tot);
  //  fprintf(stderr,"is_disp:%d, is_conductive%d, params.is_disp_ml: %d\n",is_disp,is_conductive,params.is_disp_ml);
  //if we have dispersive materials we need to create additional field variables
  auto E_nm1 = ElectricSplitField(I_tot, J_tot, K_tot);
  auto J_nm1 = CurrentDensitySplitField(I_tot, J_tot, K_tot);

  if (is_disp || inputs.params.is_disp_ml) {
    E_nm1.allocate_and_zero();
    J_nm1.allocate_and_zero();
    J_s.allocate_and_zero();
  }
  //fprintf(stderr,"Pre 14\n");
  auto J_c = CurrentDensitySplitField(I_tot, J_tot, K_tot);
  if (is_conductive) { J_c.allocate_and_zero(); }
  /*end dispersive*/

  outputs.setup_fieldsample(in_matrices["fieldsample"]);
  outputs.setup_vertex_phasors(in_matrices["campssample"], inputs.f_ex_vec.size());

  /*end of setup the output array for the sampled field*/

  /*set up the parameters for the phasor convergence procedure*/
  /*First we set dt so that an integer number of time periods fits within a sinusoidal period
   */
  double Nsteps_tmp = 0.0;
  double dt_old;
  if (inputs.params.source_mode == SourceMode::steadystate) {
    dt_old = inputs.params.dt;
    Nsteps_tmp = ceil(2. * DCPI / inputs.params.omega_an / inputs.params.dt * 3);
    inputs.params.dt = 2. * DCPI / inputs.params.omega_an * 3 / Nsteps_tmp;
  }

  //fprintf(stderr,"Pre 16\n");
  if (inputs.params.source_mode == SourceMode::steadystate && inputs.params.run_mode == RunMode::complete) {
    spdlog::info("Changed dt to {0:.10e} (was {1:.10e})", inputs.params.dt, dt_old);
  }
  Nsteps = (int) lround(Nsteps_tmp);
  //fprintf(stderr,"Pre 17\n");
  //Nsteps = (int)(floor(3*2.*DCPI/(params.omega_an*params.dt)) + 1.);//the number of time steps in a sinusoidal period
  dft_counter = 0;
  //fprintf(stderr,"Pre 18\n");
  /*params.Nt should be an integer number of Nsteps in the case of steady-state operation*/
  if (inputs.params.source_mode == SourceMode::steadystate && inputs.params.run_mode == RunMode::complete)
    if (inputs.params.Nt / Nsteps * Nsteps != inputs.params.Nt) {
      int old_Nt = inputs.params.Nt;//< For logging purposes, holds the Nt value that had to be changed
      inputs.params.Nt = inputs.params.Nt / Nsteps * Nsteps;
      spdlog::info("Changing the value of Nt to {0:d} (was {1:d})", inputs.params.Nt, old_Nt);
    }
  //fprintf(stderr,"Pre 19\n");

  if ((inputs.params.run_mode == RunMode::complete) && (inputs.params.source_mode == SourceMode::steadystate)) {
    spdlog::info("Nsteps: {0:d}", Nsteps);
  }

  /*An optimization step in the 2D (J_tot==0) case, try to work out if we have either
    of TE or TM, ie, not both*/
  int ksource_nz[4];
  for (int icomp = 0; icomp < 4; icomp++) ksource_nz[icomp] = 0;
  //fprintf(stderr,"Pre 20\n");
  if (J_tot == 0) {


    for (int icomp = 0; icomp < 4; icomp++)
      for (i = 0; i < (I_tot + 1); i++) {
        ksource_nz[icomp] = ksource_nz[icomp] ||
                            (fabs(inputs.Ksource.imag[0][i - (inputs.I0.index)][icomp]) > 1.0e-15) ||
                            (fabs(inputs.Ksource.real[0][i - (inputs.I0.index)][icomp]) > 1.0e-15);
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
    if (ksource_nz[2] || ksource_nz[1] || inputs.params.eyi_present) J_tot_bound = 1;
    else
      J_tot_bound = 0;

    //TM case
    if (ksource_nz[3] || ksource_nz[0] || inputs.params.exi_present) J_tot_p1_bound = 1;
    else
      J_tot_p1_bound = 0;
  }
  //fprintf(stderr,"Pre 22a\n");
  //fprintf(stderr,"Pre 23\n");

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
  //fprintf(stderr,"Pre 24\n");
  double time_E;
  double time_H;
  t0 = (double) time(NULL);
  //fprintf(stderr,"params.start_tind: %d\n",params.start_tind);
  //  fprintf(stdout,"params.delta.dz: %e, c: %e, params.delta.dz/c: %e\n",params.delta.dz,LIGHT_V,params.delta.dz/LIGHT_V);

  spdlog::debug("Starting main loop");

  if (TIME_MAIN_LOOP) { timers.start_timer(TimersTrackingLoop::MAIN); }

  for (unsigned int tind = inputs.params.start_tind; tind < inputs.params.Nt; tind++) {
    //fprintf(stderr,"Pos 00:\n");
    time_E = ((double) (tind + 1)) * inputs.params.dt;
    time_H = time_E - inputs.params.dt / 2.;
    //Extract phasors
    if ((dft_counter == Nsteps) && (inputs.params.run_mode == RunMode::complete)
        && (inputs.params.source_mode == SourceMode::steadystate) && inputs.params.exphasorsvolume) {

      dft_counter = 0;

      double tol = outputs.E.normalised_difference(E_copy);
      if (tol < TOL) break; //required accuracy obtained

      spdlog::debug("Phasor convergence: {} (actual) > {} (required)", tol, TOL);
      E_copy.set_values_from(outputs.E);

      outputs.E.zero();
      outputs.H.zero();
      spdlog::debug("Zeroed the phasors");

      if (inputs.params.exphasorssurface) {
        outputs.surface_phasors.zero_surface_EH();
        spdlog::debug("Zeroed the surface components");
      }
    }

    if ((inputs.params.source_mode == SourceMode::steadystate) && (inputs.params.run_mode == RunMode::complete) && inputs.params.exphasorsvolume) {

      outputs.E.set_phasors(inputs.E_s, dft_counter - 1, inputs.params.omega_an, inputs.params.dt, Nsteps);
      outputs.H.set_phasors(inputs.H_s, dft_counter, inputs.params.omega_an, inputs.params.dt, Nsteps);

      if (inputs.params.exphasorssurface) {
        for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
          outputs.surface_phasors.extractPhasorsSurface(ifx, inputs.E_s, inputs.H_s, dft_counter,
                                                        inputs.f_ex_vec[ifx] * 2 * DCPI, Nsteps, inputs.params,
                                                        inputs.params.intphasorssurface);
        }
        dft_counter++;
      }

    } else if ((inputs.params.source_mode == SourceMode::pulsed) && (inputs.params.run_mode == RunMode::complete) && inputs.params.exphasorsvolume) {
      if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        outputs.E.set_phasors(inputs.E_s, tind - 1, inputs.params.omega_an, inputs.params.dt, inputs.params.Npe);
        outputs.H.set_phasors(inputs.H_s, tind, inputs.params.omega_an, inputs.params.dt, inputs.params.Npe);
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
    if (inputs.params.source_mode == SourceMode::pulsed && inputs.params.run_mode == RunMode::complete && inputs.params.exphasorssurface) {
      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
          outputs.surface_phasors.extractPhasorsSurface(ifx, inputs.E_s, inputs.H_s, tind,
                                                        inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params.Npe,
                                                        inputs.params, inputs.params.intphasorssurface);
        }
      }
    }

    if (inputs.params.source_mode == SourceMode::pulsed && inputs.params.run_mode == RunMode::complete &&
        (outputs.vertex_phasors.there_are_vertices_to_extract_at()) &&
        ((tind - inputs.params.start_tind) % inputs.params.Np == 0)) {
      //     fprintf(stderr,"loc 01 (%d,%d,%d)\n",tind,params.start_tind,params.Np);
          //fprintf(stderr,"loc 03\n");
          //	  fprintf(stderr,"EPV 01\n");
          for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
            outputs.vertex_phasors.extractPhasorsVertices(ifx, inputs.E_s, inputs.H_s, tind, inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params);
          }
    }

    //fprintf(stderr,"Pos 02a:\n");
    if (inputs.params.source_mode == SourceMode::pulsed && inputs.params.run_mode == RunMode::complete && inputs.params.exdetintegral) {
      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        spdlog::debug("Setting Ex_t, Ey_t");

        //First need to sum up the Ex and Ey values on a plane ready for FFT, remember that Ex_t and
        // Ey_t are in row-major format whilst Exy etc. are in column major format
        for (j = inputs.params.pml.Dyl; j < (J_tot - inputs.params.pml.Dyu); j++)
          for (i = inputs.params.pml.Dxl; i < (I_tot - inputs.params.pml.Dxu); i++) {
            int m = j - inputs.params.pml.Dyl + (i - inputs.params.pml.Dxl) * (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
            Ex_t.v[m][0] = inputs.E_s.xy[inputs.params.k_det_obs][j][i] + inputs.E_s.xz[inputs.params.k_det_obs][j][i];
            Ex_t.v[m][1] = 0.;
            Ey_t.v[m][0] = inputs.E_s.yx[inputs.params.k_det_obs][j][i] + inputs.E_s.yz[inputs.params.k_det_obs][j][i];
            Ey_t.v[m][1] = 0.;
          }

        fftw_execute(Ex_t.plan);
        fftw_execute(Ey_t.plan);

        //Iterate over each mode
        for (int im = 0; im < inputs.D_tilde.num_det_modes(); im++) {

          //Now go back to column-major
          for (j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl); j++)
            for (i = 0; i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
              int m = j + i * (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
              Ex_t.cm[j][i] = Ex_t.v[m][0] + IMAGINARY_UNIT * Ex_t.v[m][1];
              Ey_t.cm[j][i] = Ey_t.v[m][0] + IMAGINARY_UNIT * Ey_t.v[m][1];
            }

          //fprintf(stderr,"Pos 02a [3]:\n");
          //Now multiply the pupil, mostly the pupil is non-zero in only a elements
          for (j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl); j++)
            for (i = 0; i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
              Ex_t.cm[j][i] *= inputs.pupil[j][i] * inputs.D_tilde.x[j][i][im];
              Ey_t.cm[j][i] *= inputs.pupil[j][i] * inputs.D_tilde.y[j][i][im];
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
                      kprop = exp(IMAGINARY_UNIT * inputs.params.z_obs * 2. * DCPI / lambda_an_t * refind *
                                  sqrt(1. - pow(lambda_an_t * inputs.f_vec.x[i] / refind, 2.) -
                                       pow(lambda_an_t * inputs.f_vec.y[j] / refind, 2.)));
                      //fprintf(stdout,"%d %d %e %e %e %e %e %e %e\n",i,j,f_vec.x[i],f_vec.y[j],real(kprop),imag(kprop),z_obs,DCPI,lambda_an_t);
                    } else {
                      kprop = exp(IMAGINARY_UNIT * (-inputs.params.air_interface + inputs.params.z_obs) * 2. * DCPI / lambda_an_t * refind *
                                  sqrt(1. - pow(lambda_an_t * inputs.f_vec.x[i] / refind, 2.) -
                                       pow(lambda_an_t * inputs.f_vec.y[j] / refind, 2.))) *
                              exp(IMAGINARY_UNIT * inputs.params.air_interface * 2. * DCPI / lambda_an_t *
                                  sqrt(1. - pow(lambda_an_t * inputs.f_vec.x[i], 2.) -
                                       pow(lambda_an_t * inputs.f_vec.y[j], 2.)));
                    }
                  } else
                    kprop = 0.;

                  Idxt += Ex_t.cm[j][i] * kprop;
                  Idyt += Ey_t.cm[j][i] * kprop;
                }
              phaseTermE = fmod(inputs.f_ex_vec[ifx] * 2. * DCPI * ((double) tind) * inputs.params.dt, 2 * DCPI);
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
        extractPhasorsPlane(iwave_lEx_bs, iwave_lEy_bs, iwave_lHx_bs, iwave_lHy_bs, inputs.E_s, inputs.H_s, I_tot,
                            J_tot, inputs.K0.index + 1, tind, inputs.params.omega_an, inputs.params.dt, inputs.params.Nt);
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
    (void)n;  // n is unused in FD derivatives – this silences the compiler warning

#pragma omp parallel default(shared) private(i, j, k, n, rho, k_loc, array_ind, Ca, Cb, Cc, alpha_l,\
                                             beta_l, gamma_l, kappa_l, sigma_l, Enp1,              \
                                             Jnp1)//,ca_vec,cb_vec,eh_vec)
    {
      n = omp_get_thread_num();
      Enp1 = 0.0;
      array_ind = 0;

      if (inputs.params.dimension == THREE || inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
          //FDTD, E_s.xy
#pragma omp for
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 1; j < J_tot; j++)
              for (i = 0; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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


                Enp1 = Ca * inputs.E_s.xy[k][j][i] + Cb * (inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i] -
                                                    inputs.H_s.zy[k][j - 1][i] - inputs.H_s.zx[k][j - 1][i]);
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.xy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dy * J_c.xy[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.xy[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xy[k][j][i];

                  E_nm1.xy[k][j][i] = inputs.E_s.xy[k][j][i];
                  J_nm1.xy[k][j][i] = J_s.xy[k][j][i];
                  J_s.xy[k][j][i] = Jnp1;

                  //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,J_s.xy[k][j][i]);
                }

                if (is_conductive && rho) { J_c.xy[k][j][i] -= rho * (Enp1 + inputs.E_s.xy[k][j][i]); }

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
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.xy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dy * J_c.xy[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.xy[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xy[k][j][i];

                  E_nm1.xy[k][j][i] = inputs.E_s.xy[k][j][i];
                  J_nm1.xy[k][j][i] = J_s.xy[k][j][i];
                  J_s.xy[k][j][i] = Jnp1;

                  //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,J_s.xy[k][j][i]);
                }

                if (is_conductive && rho) { J_c.xy[k][j][i] -= rho * (Enp1 + inputs.E_s.xy[k][j][i]); }

                eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;
                PSTD.ca[n][j - 1] = Ca;
                PSTD.cb[n][j - 1] = Cb;
              }
              if (J_tot > 1) {
                j = 0;
                eh_vec[n][j][0] = inputs.H_s.zy[k][j][i] + inputs.H_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;
                first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ey, PSTD.N_ey, inputs.E_s.xy.plan_f[n],
                                 inputs.E_s.xy.plan_b[n]);


                //fprintf(stdout,"(%d,%d) %d (of %d)\n",i,k,n,omp_get_num_threads());

                for (j = 1; j < J_tot; j++) {
                  inputs.E_s.xy[k][j][i] = PSTD.ca[n][j - 1] * inputs.E_s.xy[k][j][i] +
                                    PSTD.cb[n][j - 1] * eh_vec[n][j][0] / ((double) PSTD.N_ey);
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                Enp1 = Ca * inputs.E_s.xz[k][j][i] + Cb * (inputs.H_s.yx[k - 1][j][i] + inputs.H_s.yz[k - 1][j][i] -
                                                    inputs.H_s.yx[k][j][i] - inputs.H_s.yz[k][j][i]);
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.xz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dz * J_c.xz[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.xz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xz[k][j][i];
                  E_nm1.xz[k][j][i] = inputs.E_s.xz[k][j][i];
                  J_nm1.xz[k][j][i] = J_s.xz[k][j][i];
                  J_s.xz[k][j][i] = Jnp1;
                }

                if (is_conductive && rho) { J_c.xz[k][j][i] -= rho * (Enp1 + inputs.E_s.xz[k][j][i]); }

                inputs.E_s.xz[k][j][i] = Enp1;
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.xz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dz * J_c.xz[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.xz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.xz[k][j][i];
                  E_nm1.xz[k][j][i] = inputs.E_s.xz[k][j][i];
                  J_nm1.xz[k][j][i] = J_s.xz[k][j][i];
                  J_s.xz[k][j][i] = Jnp1;
                }

                if (is_conductive && rho) { J_c.xz[k][j][i] -= rho * (Enp1 + inputs.E_s.xz[k][j][i]); }

                eh_vec[n][k][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
                eh_vec[n][k][1] = 0.;
                PSTD.ca[n][k - 1] = Ca;
                PSTD.cb[n][k - 1] = Cb;
              }
              k = 0;
              eh_vec[n][k][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
              eh_vec[n][k][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ez, PSTD.N_ez, inputs.E_s.xz.plan_f[n],
                               inputs.E_s.xz.plan_b[n]);

              for (k = 1; k < K_tot; k++) {
                inputs.E_s.xz[k][j][i] = PSTD.ca[n][k - 1] * inputs.E_s.xz[k][j][i] -
                                  PSTD.cb[n][k - 1] * eh_vec[n][k][0] / ((double) PSTD.N_ez);
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
                if (inputs.params.is_structure) {
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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


                Enp1 = Ca * inputs.E_s.yx[k][j][i] + Cb * (inputs.H_s.zx[k][j][i - 1] + inputs.H_s.zy[k][j][i - 1] -
                                                    inputs.H_s.zx[k][j][i] - inputs.H_s.zy[k][j][i]);
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.yx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dx * J_c.yx[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.yx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx[k][j][i];
                  E_nm1.yx[k][j][i] = inputs.E_s.yx[k][j][i];
                  J_nm1.yx[k][j][i] = J_s.yx[k][j][i];
                  J_s.yx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yx[k][j][i] -= rho * (Enp1 + inputs.E_s.yx[k][j][i]); }

                inputs.E_s.yx[k][j][i] = Enp1;
              }
          //FDTD, E_s.yx
        } else {
#pragma omp for
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot_bound; j++) {
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure) {
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.yx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dx * J_c.yx[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.yx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yx[k][j][i];
                  E_nm1.yx[k][j][i] = inputs.E_s.yx[k][j][i];
                  J_nm1.yx[k][j][i] = J_s.yx[k][j][i];
                  J_s.yx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yx[k][j][i] -= rho * (Enp1 + inputs.E_s.yx[k][j][i]); }

                eh_vec[n][i][0] = inputs.H_s.zx[k][j][i] + inputs.H_s.zy[k][j][i];
                eh_vec[n][i][1] = 0.;
                PSTD.ca[n][i - 1] = Ca;
                PSTD.cb[n][i - 1] = Cb;
              }
              i = 0;
              eh_vec[n][i][0] = inputs.H_s.zx[k][j][i] + inputs.H_s.zy[k][j][i];
              eh_vec[n][i][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ex, PSTD.N_ex, inputs.E_s.yx.plan_f[n],
                               inputs.E_s.yx.plan_b[n]);

              for (i = 1; i < I_tot; i++) {
                inputs.E_s.yx[k][j][i] = PSTD.ca[n][i - 1] * inputs.E_s.yx[k][j][i] -
                                  PSTD.cb[n][i - 1] * eh_vec[n][i][0] / ((double) PSTD.N_ex);
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                Enp1 = Ca * inputs.E_s.yz[k][j][i] + Cb * (inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i] -
                                                    inputs.H_s.xy[k - 1][j][i] - inputs.H_s.xz[k - 1][j][i]);
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.yz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dz * J_c.yz[k][j][i];

                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.yz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz[k][j][i];
                  E_nm1.yz[k][j][i] = inputs.E_s.yz[k][j][i];
                  J_nm1.yz[k][j][i] = J_s.yz[k][j][i];
                  J_s.yz[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yz[k][j][i] -= rho * (Enp1 + inputs.E_s.yz[k][j][i]); }

                inputs.E_s.yz[k][j][i] = Enp1;
              }
          //FDTD, E_s.yz
        } else {
#pragma omp for
          for (j = 0; j < J_tot_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (k = 1; k < K_tot; k++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.z[k_loc];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.yz[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dz *
                                  ((1 + alpha_l) * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dz * J_c.yz[k][j][i];

                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.yz[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.yz[k][j][i];
                  E_nm1.yz[k][j][i] = inputs.E_s.yz[k][j][i];
                  J_nm1.yz[k][j][i] = J_s.yz[k][j][i];
                  J_s.yz[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.yz[k][j][i] -= rho * (Enp1 + inputs.E_s.yz[k][j][i]); }

                eh_vec[n][k][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
                eh_vec[n][k][1] = 0.;
                PSTD.ca[n][k - 1] = Ca;
                PSTD.cb[n][k - 1] = Cb;
              }
              k = 0;
              eh_vec[n][k][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
              eh_vec[n][k][1] = 0.;
              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ez, PSTD.N_ez, inputs.E_s.yz.plan_f[n],
                               inputs.E_s.yz.plan_b[n]);


              for (k = 1; k < K_tot; k++) {
                inputs.E_s.yz[k][j][i] = PSTD.ca[n][k - 1] * inputs.E_s.yz[k][j][i] +
                                  PSTD.cb[n][k - 1] * eh_vec[n][k][0] / ((double) PSTD.N_ez);
                //E_s.yz[k][j][i] = Enp1;
              }
            }
          //PSTD, E_s.yz
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
      }  //if(params.dimension==THREE || params.dimension==TE)

      //fprintf(stderr,"Pos 07:\n");
      if (inputs.params.dimension == THREE || inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
#pragma omp for
          //E_s.zx updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++)
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                Enp1 = Ca * inputs.E_s.zx[k][j][i] + Cb * (inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i] -
                                                    inputs.H_s.yx[k][j][i - 1] - inputs.H_s.yz[k][j][i - 1]);
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.zx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dx * J_c.zx[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
                  E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
                  J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                  J_s.zx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]); }

                inputs.E_s.zx[k][j][i] = Enp1;
              }
          //FDTD, E_s.zx
        } else {
#pragma omp for
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_p1_bound; j++) {
              for (i = 1; i < I_tot; i++) {
                rho = 0.;
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.x[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.zx[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dx *
                                  ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dx * J_c.zx[k][j][i];
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
                  E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
                  J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                  J_s.zx[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]); }

                eh_vec[n][i][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
                eh_vec[n][i][1] = 0.;
                PSTD.ca[n][i - 1] = Ca;
                PSTD.cb[n][i - 1] = Cb;
              }
              i = 0;
              eh_vec[n][i][0] = inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i];
              eh_vec[n][i][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ex, PSTD.N_ex, inputs.E_s.zx.plan_f[n],
                               inputs.E_s.zx.plan_b[n]);

              for (i = 1; i < I_tot; i++) {
                inputs.E_s.zx[k][j][i] = PSTD.ca[n][i - 1] * inputs.E_s.zx[k][j][i] +
                                  PSTD.cb[n][i - 1] * eh_vec[n][i][0] / ((double) PSTD.N_ex);
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
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                  if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) && (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + K - 1;
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
                if (is_conductive) rho = inputs.rho_cond.x[i];
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


              if (is_disp || inputs.params.is_disp_ml) {
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

              Enp1 = Ca * inputs.E_s.zx[k][j][i] + Cb * (inputs.H_s.yx[k][j][i] + inputs.H_s.yz[k][j][i] -
                                                  inputs.H_s.yx[k][j][i - 1] - inputs.H_s.yz[k][j][i - 1]);
              if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zx[k][j][i] -
                        1. / 2. * Cb * inputs.params.delta.dx *
                                ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
              if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dx * J_c.zx[k][j][i];

              if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                       kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zx[k][j][i];
                E_nm1.zx[k][j][i] = inputs.E_s.zx[k][j][i];
                J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                J_s.zx[k][j][i] = Jnp1;
              }
              if (is_conductive && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + inputs.E_s.zx[k][j][i]); }

              inputs.E_s.zx[k][j][i] = Enp1;
            }
      }
      //fprintf(stderr,"Pos 08:\n");
      if (inputs.params.dimension == THREE || inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
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
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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


                Enp1 = Ca * inputs.E_s.zy[k][j][i] + Cb * (inputs.H_s.xy[k][j - 1][i] + inputs.H_s.xz[k][j - 1][i] -
                                                    inputs.H_s.xy[k][j][i] - inputs.H_s.xz[k][j][i]);
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.zy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dy * J_c.zy[k][j][i];

                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
                  E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
                  J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                  J_s.zy[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]); }
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
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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
                  if (is_conductive) rho = inputs.rho_cond.y[array_ind];
                }

                alpha_l = 0.;
                beta_l = 0.;
                gamma_l = 0.;
                kappa_l = 1.;
                sigma_l = 0.;

                if (is_disp || inputs.params.is_disp_ml) {
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
                if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                  Enp1 += Cc * E_nm1.zy[k][j][i] -
                          1. / 2. * Cb * inputs.params.delta.dy *
                                  ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
                if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dy * J_c.zy[k][j][i];

                if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                  Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                         kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                  Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
                  E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
                  J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                  J_s.zy[k][j][i] = Jnp1;
                }
                if (is_conductive && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]); }

                eh_vec[n][j][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;
                PSTD.ca[n][j - 1] = Ca;
                PSTD.cb[n][j - 1] = Cb;
              }
              if (J_tot > 1) {
                j = 0;
                eh_vec[n][j][0] = inputs.H_s.xy[k][j][i] + inputs.H_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;
                first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_ey, PSTD.N_ey, inputs.E_s.zy.plan_f[n],
                                 inputs.E_s.zy.plan_b[n]);
              }
              for (j = 1; j < J_tot; j++) {
                inputs.E_s.zy[k][j][i] = PSTD.ca[n][j - 1] * inputs.E_s.zy[k][j][i] -
                                  PSTD.cb[n][j - 1] * eh_vec[n][j][0] / ((double) PSTD.N_ey);
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
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                  if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) && (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + K - 1;
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
                if (is_conductive) rho = inputs.rho_cond.y[array_ind];
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

              if (is_disp || inputs.params.is_disp_ml) {
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


              Enp1 = Ca * inputs.E_s.zy[k][j][i] + Cb * (inputs.H_s.xy[k][j - 1][i] + inputs.H_s.xz[k][j - 1][i] -
                                                  inputs.H_s.xy[k][j][i] - inputs.H_s.xz[k][j][i]);
              if ((is_disp || inputs.params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zy[k][j][i] -
                        1. / 2. * Cb * inputs.params.delta.dy *
                                ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
              if (is_conductive && rho) Enp1 += Cb * inputs.params.delta.dy * J_c.zy[k][j][i];

              if ((is_disp || inputs.params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                       kappa_l * gamma_l / (2. * inputs.params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                Jnp1 += sigma_l / EPSILON0 * gamma_l * inputs.E_s.zy[k][j][i];
                E_nm1.zy[k][j][i] = inputs.E_s.zy[k][j][i];
                J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                J_s.zy[k][j][i] = Jnp1;
              }
              if (is_conductive && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + inputs.E_s.zy[k][j][i]); }

              inputs.E_s.zy[k][j][i] = Enp1;
            }
      }
    }//end of parallel section
    //fprintf(stderr,"Pos 09:\n");
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
    /********************/

    //update terms for self consistency across scattered/total interface - E updates##
    if (inputs.params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase = exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_H, 2. * DCPI));
      double commonAmplitude = linearRamp(time_H, 1. / (inputs.params.omega_an / (2 * DCPI)), ramp_width);
      for (k = (inputs.K0.index); k <= (inputs.K1.index); k++)
        for (j = (inputs.J0.index); j <= (inputs.J1.index); j++) {
          if (inputs.I0.apply) {//Perform across I0

            if (!inputs.params.is_multilayer) array_ind = inputs.I0.index;
            else
              array_ind = (I_tot + 1) * k + inputs.I0.index;

            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
              inputs.E_s.zx[k][j][inputs.I0.index] =
                      inputs.E_s.zx[k][j][inputs.I0.index] -
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][2] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][2]));
              if (is_conductive)
                J_c.zx[k][j][inputs.I0.index] +=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][2] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][2]));
              if (inputs.params.is_disp_ml)
                J_s.zx[k][j][inputs.I0.index] +=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][2] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][2]));
            }
            if (j < (inputs.J1.index)) {
              inputs.E_s.yx[k][j][inputs.I0.index] =
                      inputs.E_s.yx[k][j][inputs.I0.index] +
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][3] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][3]));
              if (is_conductive)
                J_c.yx[k][j][inputs.I0.index] -=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][3] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][3]));
              if (inputs.params.is_disp_ml)
                J_s.yx[k][j][inputs.I0.index] -=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][3] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][3]));
            }
          }
          if (inputs.I1.apply) {//Perform across I1

            if (!inputs.params.is_multilayer) array_ind = inputs.I1.index;
            else
              array_ind = (I_tot + 1) * k + inputs.I1.index;

            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
              inputs.E_s.zx[k][j][inputs.I1.index] =
                      inputs.E_s.zx[k][j][inputs.I1.index] +
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][6] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][6]));
              if (is_conductive)
                J_c.zx[k][j][inputs.I1.index] -=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][6] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][6]));
              if (inputs.params.is_disp_ml)
                J_s.zx[k][j][inputs.I1.index] -=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][6] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][6]));
            }
            if (j < (inputs.J1.index)) {
              inputs.E_s.yx[k][j][inputs.I1.index] =
                      inputs.E_s.yx[k][j][inputs.I1.index] -
                      inputs.C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][7] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][7]));
              if (is_conductive)
                J_c.yx[k][j][inputs.I1.index] +=
                        inputs.rho_cond.x[array_ind] * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][7] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][7]));
              if (inputs.params.is_disp_ml)
                J_s.yx[k][j][inputs.I1.index] +=
                        inputs.matched_layer.kappa.x[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][7] +
                              IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][7]));
            }
          }
        }

      for (k = (inputs.K0.index); k <= (inputs.K1.index); k++)
        for (i = (inputs.I0.index); i <= (inputs.I1.index); i++) {
          if (inputs.J0.apply) {//Perform across J0
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {

              if (!inputs.params.is_multilayer) array_ind = inputs.J0.index;
              else
                array_ind = (J_tot + 1) * k + inputs.J0.index;

              inputs.E_s.zy[k][(inputs.J0.index)][i] =
                      inputs.E_s.zy[k][(inputs.J0.index)][i] +
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][2] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][2]));
              if (is_conductive)
                J_c.zy[k][(inputs.J0.index)][i] -=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][2]));
              if (inputs.params.is_disp_ml)
                J_s.zy[k][(inputs.J0.index)][i] -=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][2]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xy[k][(inputs.J0.index)][i] =
                      inputs.E_s.xy[k][(inputs.J0.index)][i] -
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][3] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][3]));
              if (is_conductive)
                J_c.xy[k][(inputs.J0.index)][i] +=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][3]));
              if (inputs.params.is_disp_ml)
                J_s.xy[k][(inputs.J0.index)][i] +=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][3]));
            }
          }
          if (inputs.J1.apply) {//Perform across J1

            if (!inputs.params.is_multilayer) array_ind = inputs.J1.index;
            else
              array_ind = (J_tot + 1) * k + inputs.J1.index;

            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC) {
              inputs.E_s.zy[k][(inputs.J1.index)][i] =
                      inputs.E_s.zy[k][(inputs.J1.index)][i] -
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][6] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][6]));
              if (is_conductive)
                J_c.zy[k][(inputs.J1.index)][i] +=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][6]));
              if (inputs.params.is_disp_ml)
                J_s.zy[k][(inputs.J1.index)][i] -=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][6]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xy[k][(inputs.J1.index)][i] =
                      inputs.E_s.xy[k][(inputs.J1.index)][i] +
                      inputs.C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][7] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][7]));
              if (is_conductive)
                J_c.xy[k][(inputs.J1.index)][i] -=
                        inputs.rho_cond.y[array_ind] * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][7]));
              if (inputs.params.is_disp_ml)
                J_s.xy[k][(inputs.J1.index)][i] +=
                        inputs.matched_layer.kappa.y[array_ind] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][7]));
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
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][2]));
              if (is_conductive)
                J_c.yz[(inputs.K0.index)][j][i] +=
                        inputs.rho_cond.z[(inputs.K0.index)] * inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][2]));
              if (inputs.params.is_disp_ml)
                J_s.yz[(inputs.K0.index)][j][i] -=
                        inputs.matched_layer.kappa.z[(inputs.K0.index)] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][2]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xz[(inputs.K0.index)][j][i] =
                      inputs.E_s.xz[(inputs.K0.index)][j][i] +
                      inputs.C.b.z[inputs.K0.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][3]));
              if (is_conductive)
                J_c.xz[(inputs.K0.index)][j][i] -=
                        inputs.rho_cond.z[(inputs.K0.index)] * inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][3]));
              if (inputs.params.is_disp_ml)
                J_s.xz[(inputs.K0.index)][j][i] +=
                        inputs.matched_layer.kappa.z[(inputs.K0.index)] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.z[inputs.K0.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][3]));
            }
          }
          if (inputs.K1.apply) {//Perform across K1
            if (j < (inputs.J1.index)) {
              inputs.E_s.yz[(inputs.K1.index)][j][i] =
                      inputs.E_s.yz[(inputs.K1.index)][j][i] +
                      inputs.C.b.z[inputs.K1.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][6] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][6]));
              if (is_conductive)
                J_c.yz[(inputs.K1.index)][j][i] -=
                        inputs.rho_cond.z[(inputs.K1.index)] * inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][6]));
              if (inputs.params.is_disp_ml)
                J_s.yz[(inputs.K1.index)][j][i] +=
                        inputs.matched_layer.kappa.z[(inputs.K1.index)] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][6] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][6]));
            }
            if (i < (inputs.I1.index)) {
              inputs.E_s.xz[(inputs.K1.index)][j][i] =
                      inputs.E_s.xz[(inputs.K1.index)][j][i] -
                      inputs.C.b.z[inputs.K1.index] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][7] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][7]));
              if (is_conductive)
                J_c.xz[(inputs.K1.index)][j][i] +=
                        inputs.rho_cond.z[(inputs.K1.index)] * inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][7]));
              if (inputs.params.is_disp_ml)
                J_s.xz[(inputs.K1.index)][j][i] -=
                        inputs.matched_layer.kappa.z[(inputs.K1.index)] * inputs.matched_layer.gamma[k] / (2. * inputs.params.dt) * inputs.C.b.z[inputs.K1.index] *
                        real(commonAmplitude * commonPhase *
                             (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][7] +
                              IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][7]));
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
                               exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
                          exp(-1.0 * DCPI *
                              pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
          //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
          if (is_conductive)
            J_c.yz[inputs.K0.index][j][i] +=
                    inputs.rho_cond.z[inputs.K0.index] * inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[0][i - (inputs.I0.index)][2] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][2]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
                    exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
          //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
          if (inputs.params.is_disp_ml) {
            J_s.yz[inputs.K0.index][j][i] -=
                    inputs.matched_layer.kappa.z[inputs.K0.index] * inputs.matched_layer.gamma[inputs.K0.index] / (2. * inputs.params.dt) *
                    inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[0][i - (inputs.I0.index)][2] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][2]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
                    exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
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
                            real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                                  IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][2]) *
                                 (-1.0 * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
                            exp(-1.0 * DCPI *
                                pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
            //E_s.yz[(int)K0[0]][j][i] = E_s.yz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (is_conductive)
              J_c.yz[inputs.K0.index][j][i] +=
                      inputs.rho_cond.z[inputs.K0.index] * inputs.C.b.z[inputs.K0.index] *
                      real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                            IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][2]) *
                           (-1.0 * IMAGINARY_UNIT) *
                           exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
                      exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
            //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            if (inputs.params.is_disp_ml) {
              J_s.yz[inputs.K0.index][j][i] -=
                      inputs.matched_layer.kappa.z[inputs.K0.index] * inputs.matched_layer.gamma[inputs.K0.index] / (2. * inputs.params.dt) *
                      inputs.C.b.z[inputs.K0.index] *
                      real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][2] +
                            IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][2]) *
                           (-1.0 * IMAGINARY_UNIT) *
                           exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
                      exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
              //J_s.yz[(int)K0[0]][j][i] -= matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
            }
          }
      for (j = 0; j < (J_tot + 1); j++)
        for (i = 0; i < I_tot; i++) {
          inputs.E_s.xz[inputs.K0.index][j][i] =
                  inputs.E_s.xz[inputs.K0.index][j][i] +
                  inputs.C.b.z[inputs.K0.index] *
                          real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                                IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][3]) *
                               (-1.0 * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2 * DCPI))) *
                          exp(-1.0 * DCPI *
                              pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
          //E_s.xz[(int)K0[0]][j][i] = E_s.xz[(int)K0[0]][j][i] + C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
          if (is_conductive)
            J_c.xz[inputs.K0.index][j][i] -=
                    inputs.rho_cond.z[inputs.K0.index] * inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][3]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2 * DCPI))) *
                    exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
          //J_c.xz[(int)K0[0]][j][i] -= rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
          if (inputs.params.is_disp_ml)
            J_s.xz[inputs.K0.index][j][i] +=
                    inputs.matched_layer.kappa.z[inputs.K0.index] * inputs.matched_layer.gamma[inputs.K0.index] / (2. * inputs.params.dt) *
                    inputs.C.b.z[inputs.K0.index] *
                    real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][3] +
                          IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][3]) *
                         (-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2 * DCPI))) *
                    exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
          //J_s.xz[(int)K0[0]][j][i] += matched_layer.kappa.z[(int)K0[0]]*matched_layer.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + IMAGINARY_UNIT*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2 ));
        }
      //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
      outputs.H.ft = real((-1.0 * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_H - inputs.params.to_l), 2. * DCPI))) *
             exp(-1.0 * DCPI * pow((time_H - inputs.params.to_l + inputs.params.delta.dz / LIGHT_V / 2.) / (inputs.params.hwhm), 2));
      //fth = real((-1.0*IMAGINARY_UNIT)*exp(-IMAGINARY_UNIT*fmod(params.omega_an*(time_H - params.to_l),2.*DCPI)))*exp( -1.0*DCPI*pow((time_H - params.to_l)/(params.hwhm),2));
    }
    //fprintf(stderr,"Pos 10:\n");

    //end of source terms
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    /********************/
    //begin parallel
#pragma omp parallel default(shared) private(i, j, k, n, k_loc,                                    \
                                             array_ind)//,ca_vec,cb_vec,eh_vec)
    {
      n = omp_get_thread_num();

      if (inputs.params.dimension == THREE || inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
//FDTD, H_s.xz
#pragma omp for
          //H_s.xz updates
          for (k = 0; k < K_tot; k++)
            for (j = 0; j < J_tot_bound; j++)
              for (i = 0; i < (I_tot + 1); i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }

                if (!inputs.materials[k][j][i])
                  inputs.H_s.xz[k][j][i] = inputs.D.a.z[k_loc] * inputs.H_s.xz[k][j][i] +
                                    inputs.D.b.z[k_loc] * (inputs.E_s.yx[k + 1][j][i] + inputs.E_s.yz[k + 1][j][i] -
                                                    inputs.E_s.yx[k][j][i] - inputs.E_s.yz[k][j][i]);
                else
                  inputs.H_s.xz[k][j][i] = inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1] * inputs.H_s.xz[k][j][i] +
                                    inputs.Dmaterial.b.z[inputs.materials[k][j][i] - 1] *
                                            (inputs.E_s.yx[k + 1][j][i] + inputs.E_s.yz[k + 1][j][i] -
                                             inputs.E_s.yx[k][j][i] - inputs.E_s.yz[k][j][i]);
              }
          //FDTD, H_s.xz
        } else {
#pragma omp for
          //H_s.xz updates
          for (j = 0; j < J_tot_bound; j++)
            for (i = 0; i < (I_tot + 1); i++) {
              for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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

                eh_vec[n][k][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
                eh_vec[n][k][1] = 0.;
              }
              k = K_tot;
              eh_vec[n][k][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
              eh_vec[n][k][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hz, PSTD.N_hz, inputs.H_s.xz.plan_f[n],
                               inputs.H_s.xz.plan_b[n]);

              for (k = 0; k < K_tot; k++) {
                inputs.H_s.xz[k][j][i] = PSTD.ca[n][k] * inputs.H_s.xz[k][j][i] +
                                  PSTD.cb[n][k] * eh_vec[n][k][0] / ((double) PSTD.N_hz);
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
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.xy[k][j][i] = inputs.D.a.y[array_ind] * inputs.H_s.xy[k][j][i] +
                                    inputs.D.b.y[array_ind] * (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] -
                                                        inputs.E_s.zy[k][j + 1][i] - inputs.E_s.zx[k][j + 1][i]);
                else
                  inputs.H_s.xy[k][j][i] = inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] * inputs.H_s.xy[k][j][i] +
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
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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

                eh_vec[n][j][0] = inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i];
                eh_vec[n][j][1] = 0.;
              }
              j = J_tot;
              eh_vec[n][j][0] = inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i];
              eh_vec[n][j][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hy, PSTD.N_hy, inputs.H_s.xy.plan_f[n],
                               inputs.H_s.xy.plan_b[n]);

              for (j = 0; j < J_tot; j++) {
                inputs.H_s.xy[k][j][i] = PSTD.ca[n][j] * inputs.H_s.xy[k][j][i] -
                                  PSTD.cb[n][j] * eh_vec[n][j][0] / ((double) PSTD.N_hy);
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
            for (j = 0; j < J_tot_p1_bound; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.yx[k][j][i] = inputs.D.a.x[array_ind] * inputs.H_s.yx[k][j][i] +
                                    inputs.D.b.x[array_ind] * (inputs.E_s.zx[k][j][i + 1] + inputs.E_s.zy[k][j][i + 1] -
                                                        inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j][i]);
                else {
                  inputs.H_s.yx[k][j][i] = inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] * inputs.H_s.yx[k][j][i] +
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
            for (j = 0; j < J_tot_p1_bound; j++) {
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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

                eh_vec[n][i][0] = inputs.E_s.zx[k][j][i] + inputs.E_s.zy[k][j][i];
                eh_vec[n][i][1] = 0.;
              }
              i = I_tot;
              eh_vec[n][i][0] = inputs.E_s.zx[k][j][i] + inputs.E_s.zy[k][j][i];
              eh_vec[n][i][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hx, PSTD.N_hx, inputs.H_s.yx.plan_f[n],
                               inputs.H_s.yx.plan_b[n]);

              for (i = 0; i < I_tot; i++) {
                inputs.H_s.yx[k][j][i] = PSTD.ca[n][i] * inputs.H_s.yx[k][j][i] +
                                  PSTD.cb[n][i] * eh_vec[n][i][0] / ((double) PSTD.N_hx);
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.materials[k][j][i]) {
                  /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                  inputs.H_s.yz[k][j][i] = inputs.D.a.z[k_loc] * inputs.H_s.yz[k][j][i] +
                                    inputs.D.b.z[k_loc] * (inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i] -
                                                    inputs.E_s.xy[k + 1][j][i] - inputs.E_s.xz[k + 1][j][i]);
                } else {
                  /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                  inputs.H_s.yz[k][j][i] = inputs.Dmaterial.a.z[inputs.materials[k][j][i] - 1] * inputs.H_s.yz[k][j][i] +
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
          for (j = 0; j < J_tot_p1_bound; j++)
#pragma omp for
            for (i = 0; i < I_tot; i++) {
              for (k = 0; k < K_tot; k++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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

                eh_vec[n][k][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
                eh_vec[n][k][1] = 0.;
              }
              k = K_tot;
              eh_vec[n][k][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
              eh_vec[n][k][1] = 0.;

              /*
    if( i==12 & j==12 ){
    for(k=0;k<K_tot;k++)
    fprintf(stdout,"%.10e ",eh_vec[n][k][0]);
    fprintf(stdout,"\n");
    }
        */

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hz, PSTD.N_hz, inputs.H_s.yz.plan_f[n],
                               inputs.H_s.yz.plan_b[n]);

              for (k = 0; k < K_tot; k++) {
                inputs.H_s.yz[k][j][i] = PSTD.ca[n][k] * inputs.H_s.yz[k][j][i] -
                                  PSTD.cb[n][k] * eh_vec[n][k][0] / ((double) PSTD.N_hz);
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
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                  if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) && (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = j;
              else
                array_ind = (J_tot + 1) * k_loc + j;
              if (!inputs.materials[k][j][i])
                inputs.H_s.xy[k][j][i] = inputs.D.a.y[array_ind] * inputs.H_s.xy[k][j][i] +
                                  inputs.D.b.y[array_ind] * (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] -
                                                    inputs.E_s.zy[k][j + 1][i] - inputs.E_s.zx[k][j + 1][i]);
              else
                inputs.H_s.xy[k][j][i] = inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] * inputs.H_s.xy[k][j][i] +
                                  inputs.Dmaterial.b.y[inputs.materials[k][j][i] - 1] *
                                          (inputs.E_s.zy[k][j][i] + inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j + 1][i] -
                                           inputs.E_s.zx[k][j + 1][i]);
            }

#pragma omp for
        //H_s.yx update
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < (J_tot + 1); j++)
            for (i = 0; i < I_tot; i++) {
              k_loc = k;
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                  if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) && (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + K - 1;
                  else
                    k_loc = inputs.params.pml.Dzl + 1;
                }
              if (!inputs.params.is_multilayer) array_ind = i;
              else
                array_ind = (I_tot + 1) * k_loc + i;
              if (!inputs.materials[k][j][i])
                inputs.H_s.yx[k][j][i] = inputs.D.a.x[array_ind] * inputs.H_s.yx[k][j][i] +
                                  inputs.D.b.x[array_ind] * (inputs.E_s.zx[k][j][i + 1] + inputs.E_s.zy[k][j][i + 1] -
                                                    inputs.E_s.zx[k][j][i] - inputs.E_s.zy[k][j][i]);
              else
                inputs.H_s.yx[k][j][i] = inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] * inputs.H_s.yx[k][j][i] +
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

      if (inputs.params.dimension == THREE || inputs.params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        if (solver_method == SolverMethod::FiniteDifference) {
//FDTD, H_s.zy
#pragma omp for
          //H_s.zy update
          for (k = 0; k < (K_tot + 1); k++)
            for (j = 0; j < J_tot; j++)
              for (i = 0; i < I_tot; i++) {
                k_loc = k;
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = j;
                else
                  array_ind = (J_tot + 1) * k_loc + j;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.zy[k][j][i] = inputs.D.a.y[array_ind] * inputs.H_s.zy[k][j][i] +
                                    inputs.D.b.y[array_ind] * (inputs.E_s.xy[k][j + 1][i] + inputs.E_s.xz[k][j + 1][i] -
                                                        inputs.E_s.xy[k][j][i] - inputs.E_s.xz[k][j][i]);
                else
                  inputs.H_s.zy[k][j][i] = inputs.Dmaterial.a.y[inputs.materials[k][j][i] - 1] * inputs.H_s.zy[k][j][i] +
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
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
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

                eh_vec[n][j][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
                eh_vec[n][j][1] = 0.;
              }
              j = J_tot;
              eh_vec[n][j][0] = inputs.E_s.xy[k][j][i] + inputs.E_s.xz[k][j][i];
              eh_vec[n][j][1] = 0.;

              first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hy, PSTD.N_hy, inputs.H_s.zy.plan_f[n],
                               inputs.H_s.zy.plan_b[n]);

              for (j = 0; j < J_tot; j++) {
                inputs.H_s.zy[k][j][i] = PSTD.ca[n][j] * inputs.H_s.zy[k][j][i] +
                                  PSTD.cb[n][j] * eh_vec[n][j][0] / ((double) PSTD.N_hy);
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
                if (inputs.params.is_structure)
                  if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                    if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) &&
                        (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                      k_loc = k - inputs.structure[i][1];
                    else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                      k_loc = inputs.params.pml.Dzl + K - 1;
                    else
                      k_loc = inputs.params.pml.Dzl + 1;
                  }
                if (!inputs.params.is_multilayer) array_ind = i;
                else
                  array_ind = (I_tot + 1) * k_loc + i;
                if (!inputs.materials[k][j][i])
                  inputs.H_s.zx[k][j][i] = inputs.D.a.x[array_ind] * inputs.H_s.zx[k][j][i] +
                                    inputs.D.b.x[array_ind] * (inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i] -
                                                        inputs.E_s.yx[k][j][i + 1] - inputs.E_s.yz[k][j][i + 1]);
                else
                  inputs.H_s.zx[k][j][i] = inputs.Dmaterial.a.x[inputs.materials[k][j][i] - 1] * inputs.H_s.zx[k][j][i] +
                                    inputs.Dmaterial.b.x[inputs.materials[k][j][i] - 1] *
                                            (inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i] -
                                             inputs.E_s.yx[k][j][i + 1] - inputs.E_s.yz[k][j][i + 1]);
              }
          //FDTD, H_s.zx
        } else {
#pragma omp for
        //H_s.zx update
        for (k = 0; k < (K_tot + 1); k++)
          for (j = 0; j < J_tot_bound; j++) {
            for (i = 0; i < I_tot; i++) {
              k_loc = k;
              if (inputs.params.is_structure)
                if (k > inputs.params.pml.Dzl && k < (inputs.params.pml.Dzl + K)) {
                  if ((k - inputs.structure[i][1]) < (K + inputs.params.pml.Dzl) && (k - inputs.structure[i][1]) > inputs.params.pml.Dzl)
                    k_loc = k - inputs.structure[i][1];
                  else if ((k - inputs.structure[i][1]) >= (K + inputs.params.pml.Dzl))
                    k_loc = inputs.params.pml.Dzl + K - 1;
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

              eh_vec[n][i][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
              eh_vec[n][i][1] = 0.;
            }
            i = I_tot;
            eh_vec[n][i][0] = inputs.E_s.yx[k][j][i] + inputs.E_s.yz[k][j][i];
            eh_vec[n][i][1] = 0.;


            first_derivative(eh_vec[n], eh_vec[n], PSTD.dk_hx,
                             PSTD.N_hx, inputs.H_s.zx.plan_f[n], inputs.H_s.zx.plan_b[n]);

            for (i = 0; i < I_tot; i++) {
              inputs.H_s.zx[k][j][i] = PSTD.ca[n][i] * inputs.H_s.zx[k][j][i] -
                                PSTD.cb[n][i] *
                                        eh_vec[n][i][0] / ((double) PSTD.N_hx);
            }
          }
//PSTD, H_s.zx
        }// if (solver_method == DerivativeMethod::FiniteDifference) (else PseudoSpectral)
      }//(params.dimension==THREE || params.dimension==TE)
    }  //end parallel
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    //fprintf(stderr,"Pos 11b:\n");
    //update terms for self consistency across scattered/total interface - E updates
    if (inputs.params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase = exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * time_E, 2. * DCPI));
      double commonAmplitude = linearRamp(time_E, 1. / (inputs.params.omega_an / (2 * DCPI)), ramp_width);
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
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][0] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][0]));
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.yx[k][j][(inputs.I0.index) - 1] =
                      inputs.H_s.yx[k][j][(inputs.I0.index) - 1] -
                      inputs.D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][1] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][1]));
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
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][4] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][4]));
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.yx[k][j][(inputs.I1.index)] =
                      inputs.H_s.yx[k][j][(inputs.I1.index)] +
                      inputs.D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Isource.real[k - (inputs.K0.index)][j - (inputs.J0.index)][5] +
                                    IMAGINARY_UNIT * inputs.Isource.imag[k - (inputs.K0.index)][j - (inputs.J0.index)][5]));
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
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][0] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][0]));

            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.xy[k][(inputs.J0.index) - 1][i] =
                      inputs.H_s.xy[k][(inputs.J0.index) - 1][i] +
                      inputs.D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][1] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][1]));
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
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][4] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][4]));
            if (k < (inputs.K1.index) || inputs.params.dimension == Dimension::TRANSVERSE_MAGNETIC)
              inputs.H_s.xy[k][(inputs.J1.index)][i] =
                      inputs.H_s.xy[k][(inputs.J1.index)][i] -
                      inputs.D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Jsource.real[k - (inputs.K0.index)][i - (inputs.I0.index)][5] +
                                    IMAGINARY_UNIT * inputs.Jsource.imag[k - (inputs.K0.index)][i - (inputs.I0.index)][5]));
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
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][0] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][0]));
            if (j < (inputs.J1.index))
              inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.xz[(inputs.K0.index) - 1][j][i] -
                      inputs.D.b.z[(inputs.K0.index) - 1] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][1] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][1]));
          }
          if (inputs.K1.apply) {//Perform across K1
            if (i < (inputs.I1.index))
              inputs.H_s.yz[(inputs.K1.index)][j][i] =
                      inputs.H_s.yz[(inputs.K1.index)][j][i] -
                      inputs.D.b.z[(inputs.K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][4] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][4]));
            if (j < (inputs.J1.index))
              inputs.H_s.xz[(inputs.K1.index)][j][i] =
                      inputs.H_s.xz[(inputs.K1.index)][j][i] +
                      inputs.D.b.z[(inputs.K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][5] +
                                    IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][5]));
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
                               exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_E - inputs.params.to_l), 2 * DCPI))) *
                          exp(-1. * DCPI * pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
          //broadband source term
          if (inputs.params.eyi_present)
            inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.xz[(inputs.K0.index) - 1][j][i] - inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.y[tind][j][i];
        }
        //fprintf(stderr,"Pos 11e\n");
        for (i = 0; i < I_tot; i++) {
          inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                  inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                  inputs.D.b.z[(inputs.K0.index) - 1] *
                          real((inputs.Ksource.real[0][i - (inputs.I0.index)][0] +
                                IMAGINARY_UNIT * inputs.Ksource.imag[0][i - (inputs.I0.index)][0]) *
                               (-1. * IMAGINARY_UNIT) *
                               exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_E - inputs.params.to_l), 2 * DCPI))) *
                          exp(-1. * DCPI * pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
          //broadband source term
          if (inputs.params.exi_present)
            inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.yz[(inputs.K0.index) - 1][j][i] + inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.x[tind][j][i];
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
                            real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][1] +
                                  IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][1]) *
                                 (-1. * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_E - inputs.params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
            //broadband source term
            if (inputs.params.eyi_present)
              inputs.H_s.xz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.xz[(inputs.K0.index) - 1][j][i] - inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.y[tind][j][i];
          }
        //fprintf(stderr,"Pos 11h\n");
        for (j = 0; j < (J_tot + 1); j++)
          for (i = 0; i < I_tot; i++) {
            inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                    inputs.H_s.yz[(inputs.K0.index) - 1][j][i] +
                    inputs.D.b.z[(inputs.K0.index) - 1] *
                            real((inputs.Ksource.real[j - (inputs.J0.index)][i - (inputs.I0.index)][0] +
                                  IMAGINARY_UNIT * inputs.Ksource.imag[j - (inputs.J0.index)][i - (inputs.I0.index)][0]) *
                                 (-1. * IMAGINARY_UNIT) *
                                 exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_E - inputs.params.to_l), 2 * DCPI))) *
                            exp(-1. * DCPI * pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
            //broadband source term
            if (inputs.params.exi_present)
              inputs.H_s.yz[(inputs.K0.index) - 1][j][i] =
                      inputs.H_s.yz[(inputs.K0.index) - 1][j][i] + inputs.D.b.z[(inputs.K0.index) - 1] * inputs.Ei.x[tind][j][i];
          }
        //fprintf(stderr,"Pos 11i\n");
      }
      outputs.E.ft = real((-1. * IMAGINARY_UNIT) * exp(-IMAGINARY_UNIT * fmod(inputs.params.omega_an * (time_E - inputs.params.to_l), 2 * DCPI))) *
             exp(-1. * DCPI * pow((time_E - inputs.params.to_l) / (inputs.params.hwhm), 2.));
      //fprintf(stderr,"Pos 11j\n");
    }
    if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }

    if (inputs.params.exphasorssurface || inputs.params.exphasorsvolume || inputs.params.exdetintegral || outputs.vertex_phasors.there_are_vertices_to_extract_at()) {
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

      if( (tind % Nsteps) == 0 ){
        outputs.E.angular_norm = 0.0;
        outputs.H.angular_norm = 0.0;

        for(int ifx=0; ifx<inputs.f_ex_vec.size(); ifx++){
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

        outputs.E.add_to_angular_norm(tind, Nsteps, inputs.params);
        outputs.H.add_to_angular_norm(tind, Nsteps, inputs.params);

        for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
          extractPhasorENorm(&E_norm[ifx], outputs.E.ft, tind, inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params.dt, Nsteps);
          extractPhasorHNorm(&H_norm[ifx], outputs.H.ft, tind, inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params.dt, Nsteps);
        }
      } else {
        if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {

          outputs.E.add_to_angular_norm(tind, inputs.params.Npe, inputs.params);
          outputs.H.add_to_angular_norm(tind, inputs.params.Npe, inputs.params);

          for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
            extractPhasorENorm(&E_norm[ifx], outputs.E.ft, tind, inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params.dt, inputs.params.Npe);
            extractPhasorHNorm(&H_norm[ifx], outputs.H.ft, tind, inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params.dt, inputs.params.Npe);
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
    if ((inputs.params.source_mode == SourceMode::steadystate) && (tind == (inputs.params.Nt - 1)) && (inputs.params.run_mode == RunMode::complete) &&
        inputs.params.exphasorsvolume) {
      fprintf(stdout, "Iteration limit reached, setting output fields to last complete DFT\n");
      outputs.E.set_values_from(E_copy);
    }
    //fprintf(stderr,"Post-iter 4\n");
    fflush(stdout);
    //fprintf(stderr,"Post-iter 5\n");
    //fprintf(stderr,"%s %d %d\n",tdfdirstr, strcmp(tdfdirstr,""),are_equal(tdfdirstr,""));
    if (inputs.params.has_tdfdir && (tind % inputs.params.Np) == 0) {
      fprintf(stderr,"Saving field\n");
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
    spdlog::info("Time elapsed in main loop: {0:e}", timers.time_ellapsed_by(TimersTrackingLoop::MAIN));
  }
  //save state of fdtdgrid

  //fprintf(stderr,"Pos 12\n");
  if (inputs.params.run_mode == RunMode::complete && inputs.params.exphasorsvolume) {
    outputs.E.normalise_volume();
    outputs.H.normalise_volume();
  }

  //fprintf(stderr,"Pos 13\n");
  if (inputs.params.run_mode == RunMode::complete && inputs.params.exphasorssurface) {
    spdlog::info("Surface phasors");
    for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
      outputs.surface_phasors.normalise_surface(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
    }
  }
  if (inputs.params.run_mode == RunMode::complete && outputs.vertex_phasors.there_are_vertices_to_extract_at()) {
    spdlog::info("Vertex phasors");
    for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
      outputs.vertex_phasors.normalise_vertices(ifx, E_norm[ifx], H_norm[ifx]);
      spdlog::info("\tE_norm[{0:d}]: {1:.5e} {2:.5e}", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
    }
  }

  //fprintf(stderr,"Pos 14\n");
  if (inputs.params.source_mode == SourceMode::pulsed && inputs.params.run_mode == RunMode::complete && inputs.params.exdetintegral) {
    for (int im = 0; im < inputs.D_tilde.num_det_modes(); im++)
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        outputs.ID.x[ifx][im] = outputs.ID.x[ifx][im] / E_norm[ifx];
        outputs.ID.y[ifx][im] = outputs.ID.y[ifx][im] / E_norm[ifx];

        outputs.ID.x_real[ifx][im] = real(outputs.ID.x[ifx][im]);
        outputs.ID.x_imag[ifx][im] = imag(outputs.ID.x[ifx][im]);

        outputs.ID.y_real[ifx][im] = real(outputs.ID.y[ifx][im]);
        outputs.ID.y_imag[ifx][im] = imag(outputs.ID.y[ifx][im]);
      }
  }

  /* Find the maximum absolute value of residual field in the grid, and write too the outputs */
  maxfield = max(inputs.E_s.largest_field_value(), inputs.H_s.largest_field_value());
  outputs.set_maxresfield(maxfield, false);

  // setup interpolated field outputs and labels (if necessary)
  outputs.setup_interpolation_outputs(inputs.params);

  /*Now export 3 matrices, a vertex list, a matrix of complex amplitudes at
    these vertices and a list of facets*/
  // if we need to extract the phasors, we will need to allocate memory in the output
  bool extracting_phasors = (inputs.params.exphasorssurface && inputs.params.run_mode == RunMode::complete);
  if (extracting_phasors) {
    //first regenerate the mesh since we threw away the facet list before iterating
    mxArray *dummy_vertex_list;
    if (J_tot == 0)
      conciseCreateBoundary(inputs.cuboid[0], inputs.cuboid[1], inputs.cuboid[4], inputs.cuboid[5], &dummy_vertex_list,
                            &mx_surface_facets);
    else
      conciseTriangulateCuboidSkip(inputs.cuboid[0], inputs.cuboid[1], inputs.cuboid[2], inputs.cuboid[3], inputs.cuboid[4], inputs.cuboid[5],
                                   inputs.params.spacing_stride, &dummy_vertex_list,
                                   &mx_surface_facets);
    mxDestroyArray(dummy_vertex_list);

    //now create and populate the vertex list
    outputs.surface_phasors.create_vertex_list(inputs.input_grid_labels);
  }
  // assign to the output
  outputs.assign_surface_phasor_outputs(!extracting_phasors, mx_surface_facets);

  /*End of FDTD iteration*/

  /* Free the additional data structures used to cast the matlab arrays*/
  if (inputs.params.exphasorssurface && inputs.params.run_mode == RunMode::complete) {
    mxDestroyArray(mx_surface_vertices);
  }

  free(E_norm);
  free(H_norm);

  if (inputs.params.source_mode == SourceMode::steadystate && inputs.params.run_mode == RunMode::complete) {
    mxDestroyArray(E_copy_MATLAB_data[0]);
    mxDestroyArray(E_copy_MATLAB_data[1]);
    mxDestroyArray(E_copy_MATLAB_data[2]);
  }
  //must destroy surface_phasors.mx_surface_amplitudes

  return outputs;
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

void extractPhasorsPlane(Matrix<complex<double>> &iwave_lEx_bs, Matrix<complex<double>> &iwave_lEy_bs,
                         Matrix<complex<double>> &iwave_lHx_bs, Matrix<complex<double>> &iwave_lHy_bs,
                         ElectricSplitField &E, MagneticSplitField &H, int I_tot, int J_tot, int K1,
                         int n, double omega, double dt, int Nt) {
  complex<double> phaseTerm = fmod(omega * ((double) n) * dt, 2 * DCPI);

  for (int j = 0; j < J_tot; j++) {
    for (int i = 0; i < (I_tot + 1); i++) {
      complex<double> subResult = 0.;
      //Eyz
      subResult = (E.yz[K1][j][i] + E.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);
      iwave_lEy_bs[i][j] += subResult;
      //Hxz
      subResult = (H.xz[K1 - 1][j][i] + H.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);
      iwave_lHx_bs[i][j] += subResult;
    }
  }
  for (int j = 0; j < (J_tot + 1); j++) {
    for (int i = 0; i < I_tot; i++) {
      complex<double> subResult = 0.;
      //Exz
      subResult = (E.xz[K1][j][i] + E.xy[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);
      iwave_lEx_bs[i][j] += subResult;
      //Hyz
      subResult = (H.yz[K1 - 1][j][i] + H.yx[K1][j][i]) * exp(phaseTerm * IMAGINARY_UNIT) * 1. / ((double) Nt);
      iwave_lHy_bs[i][j] += subResult;
    }
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
