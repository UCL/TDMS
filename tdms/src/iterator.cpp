/*****************************************************************
 *  Application.:  main time propogation algorithm
 *  Description.:  Contains the main FDTD loop as well as other functions
 *                 such as phasor extraction etc. Works in both pulsed 
 *                 and steady state mode.
 ******************************************************************/
#include <complex>
#include <algorithm>
#include <cstring>
#include <ctime>
#include <omp.h>
#include "iterator.h"
#include "mesh_base.h"
#include "numerical_derivative.h"
#include "array_init.h"
#include "globals.h"
#include "interface.h"
#include "interpolate.h"
#include "matlabio.h"
#include "shapes.h"
#include "source.h"
#include "timer.h"
#include "utils.h"


using namespace std;

//whether of not to time execution
#define TIME_EXEC false
//time main loop
#define TIME_MAIN_LOOP true
//threshold used to terminate the steady state iterations
#define TOL 1e-6
//parameter to control the with of the ramp when introducing the waveform in steady state mode
#define ramp_width 4.


/*This mex function will take in the following arguments and perform the
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

  fdtdgrid - A structre with the following members, each of which is a 3 dimensional
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

  conductive_aux - auxilliary parameters required to model conductive multilayer

  dispersive_aux - auxilliary parameters required to model dispersive multilayer

  structure - 2 x (I_tot+1) integer array describing the grating stucture, if one is present

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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  auto params = SimulationParameters();

  auto E_s = ElectricSplitField();
  auto H_s = MagneticSplitField();
  auto J_s = CurrentDensitySplitField();

  auto E = ElectricField();
  auto H = MagneticField();
  auto E_copy = ElectricField();  // Used to check convergence with E - E_copy

  double ***surface_EHr, ***surface_EHi;
  double rho;
  double alpha_l, beta_l, gamma_l;
  double kappa_l, sigma_l;
  double t0;

  double Ca, Cb, Cc;     //used by interpolation scheme
  //the C and D vars for free space and pml
  double Enp1, Jnp1;
  //these are used for boot strapping. There is currently no way of exporting this.
  double **iwave_lEx_Rbs, **iwave_lEy_Rbs, **iwave_lHx_Rbs, **iwave_lHy_Rbs, **iwave_lEx_Ibs,
          **iwave_lEy_Ibs, **iwave_lHx_Ibs, **iwave_lHy_Ibs;
  double maxfield = 0, tempfield;

  //refractive index of the first layer of the multilayer, or of the bulk of homogeneous
  double refind;

  //PSTD storage
  fftw_complex *dk_e_x, *dk_e_y, *dk_e_z, *dk_h_x, *dk_h_y, *dk_h_z;
  int N_e_x, N_e_y, N_e_z, N_h_x, N_h_y, N_h_z;
  double phaseTermE;
  complex<double> cphaseTermE;
  double lambda_an_t;
  //end PSTD storage

  uint8_t ***materials;
  double ***camplitudesR, ***camplitudesI;
  mxArray *mx_camplitudes;

  int i, j, k, n, k_loc, ndims, K;
  int input_counter = 0;
  int Nsteps = 0, dft_counter = 0;
  int **surface_vertices, n_surface_vertices = 0;
  int Ni_tdf = 0, Nk_tdf = 0;

#ifdef FDFLAG
  int skip_tdf = 6;
#else
  int skip_tdf = 1;
#endif

  const mwSize *dimptr_out;
  mwSize *dims;
  dims = (mwSize *) malloc(3 * sizeof(mwSize));
  mwSize *label_dims;
  label_dims = (mwSize *) malloc(2 * sizeof(mwSize));
  mxArray *dummy_array[3];
  mxArray *mx_surface_vertices, *mx_surface_facets, *mx_surface_amplitudes;
  mxArray *mx_Idx, *mx_Idy;
  double **Idx_re, **Idx_im, **Idy_re, **Idy_im;
  complex<double> **Idx, **Idy;
  complex<double> Idxt, Idyt, kprop;

  fprintf(stdout, "Using %d OMP threads\n", omp_get_max_threads());

  if (nrhs != 49) { throw runtime_error("Expected 49 inputs. Had " + to_string(nrhs)); }

  if (nlhs != 31) { throw runtime_error("31 outputs required. Had " + to_string(nlhs)); }

  /*Get fdtdgrid */
  assert_is_struct(prhs[input_counter], "fdtdgrid, argument " + to_string(input_counter));
  init_grid_arrays(prhs[input_counter], E_s, H_s, materials);
  int I_tot = E_s.I_tot, J_tot = E_s.J_tot, K_tot = E_s.K_tot;
  input_counter++;
  //fprintf(stderr,"Got fdtdgrid\n");

  /*Get Cmaterials */
  assert_is_struct(prhs[input_counter], "Cmaterials, argument " + to_string(input_counter));
  auto Cmaterial = CMaterial(prhs[input_counter]);
  input_counter++;
  //fprintf(stderr,"Got Cmaterials\n");

  /*Get Dmaterials */
  assert_is_struct(prhs[input_counter], "Dmaterials, argument " + to_string(input_counter));
  auto Dmaterial = DMaterial(prhs[input_counter]);
  input_counter++;
  //fprintf(stderr,"Got Dmaterials\n");

  /*Get C */
  assert_is_struct(prhs[input_counter], "C, argument " + to_string(input_counter));
  auto C = CCollection(prhs[input_counter]);
  params.is_disp_ml = C.is_disp_ml;
  params.is_multilayer = C.is_multilayer;
  input_counter++;
  //fprintf(stderr,"Got C\n");

  /*Get D */
  assert_is_struct(prhs[input_counter], "D, argument " + to_string(input_counter));
  auto D = DCollection(prhs[input_counter]);
  input_counter++;
  //fprintf(stderr,"Got D\n");

  /*Get freespace*/  // Cby Cbz Dbx Dby Dbz are unused
  assert_is_struct_with_n_fields(prhs[input_counter], 6, "freespace, argument " + to_string(input_counter));
  auto freespace_Cbx = mxGetPr(ptr_to_vector_in(prhs[input_counter], "Cbx", "freespace"));
  input_counter++;
  //fprintf(stderr,"Got freespace\n");

  /*Get disp_params */
  assert_is_struct_with_n_fields(prhs[input_counter], 3, "disp_params, argument " + to_string(input_counter));
  auto alpha = mxGetPr(ptr_to_vector_or_empty_in(prhs[input_counter], "alpha", "disp_params"));
  auto beta  = mxGetPr(ptr_to_vector_or_empty_in(prhs[input_counter], "beta",  "disp_params"));
  auto gamma = mxGetPr(ptr_to_vector_or_empty_in(prhs[input_counter], "gamma", "disp_params"));
  input_counter++;
  //fprintf(stderr,"Got disp_params\n");

  /*Get delta params*/
  assert_is_struct_with_n_fields(prhs[input_counter], 3, "delta, argument " + to_string(input_counter));
  auto dx = *mxGetPr(ptr_to_vector_in(prhs[input_counter], "x", "delta"));
  auto dy = *mxGetPr(ptr_to_vector_in(prhs[input_counter], "y", "delta"));
  auto dz = *mxGetPr(ptr_to_vector_in(prhs[input_counter], "z", "delta"));
  input_counter++;
  //fprintf(stderr,"Got delta params\n");

  /*Get interface*/
  assert_is_struct_with_n_fields(prhs[input_counter], 6, "interface, argument " + to_string(input_counter));
  auto I0 = InterfaceComponent(prhs[input_counter], "I0");
  auto I1 = InterfaceComponent(prhs[input_counter], "I1");
  auto J0 = InterfaceComponent(prhs[input_counter], "J0");
  auto J1 = InterfaceComponent(prhs[input_counter], "J1");
  auto K0 = InterfaceComponent(prhs[input_counter], "K0");
  auto K1 = InterfaceComponent(prhs[input_counter], "K1");
  input_counter++;
  //fprintf(stderr,"Got interface\n");

  auto Isource = Source(prhs[input_counter++], J1.index - J0.index + 1, K1.index - K0.index + 1, "Isource");
  auto Jsource = Source(prhs[input_counter++], I1.index - I0.index + 1, K1.index - K0.index + 1, "Jsource");
  auto Ksource = Source(prhs[input_counter++], I1.index - I0.index + 1, J1.index - J0.index + 1, "Ksource");

  /*Get grid_labels*/
  assert_is_struct_with_n_fields(prhs[input_counter], 3, "grid_labels, argument " + to_string(input_counter));
  auto input_grid_labels = GridLabels(prhs[input_counter]);
  input_counter++;
  //fprintf(stderr,"Got   grid_labels\n");

  params.omega_an = double_in(prhs[input_counter++], "omega_an");
  params.to_l = double_in(prhs[input_counter++], "to_l");
  params.hwhm = double_in(prhs[input_counter++], "hwhm");
  params.pml.Dxl = int_cast_from_double_in(prhs[input_counter++], "Dxl");
  params.pml.Dxu = int_cast_from_double_in(prhs[input_counter++], "Dxu");
  params.pml.Dyl = int_cast_from_double_in(prhs[input_counter++], "Dyl");
  params.pml.Dyu = int_cast_from_double_in(prhs[input_counter++], "Dyu");
  params.pml.Dzl = int_cast_from_double_in(prhs[input_counter++], "Dzl");
  params.pml.Dzu = int_cast_from_double_in(prhs[input_counter++], "Dzu");

  params.Nt = int_cast_from_double_in(prhs[input_counter++], "Nt");
  params.dt = double_in(prhs[input_counter++], "dt");
  params.start_tind = int_cast_from_double_in(prhs[input_counter++], "tind");

  params.set_source_mode(string_in(prhs[input_counter++], "sourcemode"));
  params.set_run_mode(string_in(prhs[input_counter++], "runmode"));

  params.exphasorsvolume = bool_cast_from_double_in(prhs[input_counter++], "exphasorsvolume");
  params.exphasorssurface = bool_cast_from_double_in(prhs[input_counter++], "exphasorssurface");
  params.intphasorssurface = bool_cast_from_double_in(prhs[input_counter++], "intphasorssurface");

  /*Get phasorsurface*/
  /*Only do if exphasorssurface is true*/
  auto cuboid = Cuboid();
  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    cuboid.initialise(prhs[input_counter], J_tot);
  }
  input_counter++;
  //fprintf(stderr,"Got   phasorsurface\n");

  params.set_phasorinc(mxGetPr(prhs[input_counter++]));
  params.set_dimension(string_in(prhs[input_counter++], "dimension"));

  /*Get conductive_aux */
  assert_is_struct_with_n_fields(prhs[input_counter], 3, "conductive_aux");
  auto rho_cond = XYZVectors();
  rho_cond.x = mxGetPr(ptr_to_vector_in(prhs[input_counter], "rho_x", "conductive_aux"));
  rho_cond.y = mxGetPr(ptr_to_vector_in(prhs[input_counter], "rho_y", "conductive_aux"));
  rho_cond.z = mxGetPr(ptr_to_vector_in(prhs[input_counter], "rho_z", "conductive_aux"));
  input_counter++;

  /*Get dispersive_aux*/
  auto ml = DispersiveMultiLayer(prhs[input_counter++]);
  
  /*Get structure*/
  auto structure = GratingStructure(prhs[input_counter++], I_tot);
  params.is_structure = structure.has_elements();

  /*Get f_ex_vec*/
  auto f_ex_vec = FrequencyExtractVector(prhs[input_counter++], params.omega_an);
  
  /*Get exdetintegral*/
  if (!mxIsEmpty(prhs[input_counter])) {
    params.exdetintegral = bool_cast_from_double_in(prhs[input_counter], "exdetintegral");
  }
  input_counter++;

  auto f_vec = FrequencyVectors();
  auto pupil = Pupil();
  auto D_tilde = DTilde();

  if (params.exdetintegral) {
    f_vec.initialise(prhs[input_counter++]);
    pupil.initialise(prhs[input_counter++], f_vec.x.size(), f_vec.y.size());
    D_tilde.initialise(prhs[input_counter++], f_vec.x.size(), f_vec.y.size());

    if (!mxIsEmpty(prhs[input_counter])) {
      params.k_det_obs = int_cast_from_double_in(prhs[input_counter], "k_det_obs") - 1;
    }
    input_counter++;

    params.z_obs = input_grid_labels.z[params.k_det_obs];
  }
  else {
    input_counter += 4; // advance beyond fields which were not read in as exdetintegral was false
  }

  /*Get air_interface*/
  if (!mxIsEmpty(prhs[input_counter])) {
    params.air_interface_present = true;
    params.air_interface = double_in(prhs[input_counter], "air_interface");
    //fprintf(stderr, "air_interface: %e\nz_obs: %e\n", params.air_interface, params.z_obs);
  }
  input_counter++;

  params.interp_mat_props = bool_cast_from_double_in(prhs[input_counter++], "intmatprops");

  /*Get intmethod*/
  if (!mxIsEmpty(prhs[input_counter])) {
    params.interp_method = InterpolationMethod(int_cast_from_double_in(prhs[input_counter], "intmethod"));
  }
  fprintf(stderr, "intmethod=%d\n", params.interp_method);
  input_counter++;

  /*Get tdfield*/
  auto Ei = IncidentField(prhs[input_counter++]);
  params.exi_present = Ei.x.has_elements();
  params.eyi_present = Ei.y.has_elements();

  /*Get tdfdir*/
  //fprintf(stderr,"tdfdir: %d (%d)\n", mxIsChar(prhs[input_counter]),input_counter);
  auto ex_td_field_exporter = TDFieldExporter2D();

  if (mxIsChar(prhs[input_counter])) {

    ex_td_field_exporter.folder_name = string_in(prhs[input_counter], "tdfdir").c_str();

    for (k = 0; k < K_tot; k++)
      if ((k % skip_tdf) == 0) Nk_tdf++;

    for (i = 0; i < I_tot; i++)
      if ((i % skip_tdf) == 0) Ni_tdf++;
    fprintf(stderr, "Ni_tdf=%d, Nk_tdf=%d\n", Ni_tdf, Nk_tdf);

    if (!are_equal(ex_td_field_exporter.folder_name, "")){
      params.has_tdfdir = true;
      ex_td_field_exporter.allocate(Ni_tdf, Nk_tdf);
    }

    input_counter++;
  }
  /*Got tdfdir*/

  auto fieldsample = FieldSample(prhs[input_counter++]);
  auto campssample = ComplexAmplitudeSample(prhs[input_counter++]);

  /*Deduce the refractive index of the first layer of the multilayer, or of the bulk of homogeneous*/
  refind = sqrt(1. / (freespace_Cbx[0] / params.dt * dx) / eo);
  fprintf(stderr, "refind=%e\n", refind);

  /*Setup temporary storage for detector sensitivity evaluation*/
  auto Ex_t = DetectorSensitivityArrays();
  auto Ey_t = DetectorSensitivityArrays();

  if (params.exdetintegral) {
    int n0 = I_tot - params.pml.Dxl - params.pml.Dxu;
    int n1 = J_tot - params.pml.Dyl - params.pml.Dyu;
    Ex_t.initialise(n1, n0);
    Ey_t.initialise(n1, n0);
  }

  CCoefficientMatrix ca_vec, cb_vec, cc_vec;
  EHVec eh_vec;

#ifndef FDFLAG// only perform if using the PSTD method
  int max_IJK = E_s.max_IJK_tot(), n_threads = omp_get_max_threads();
  ca_vec.allocate(n_threads, max_IJK + 1);
  cb_vec.allocate(n_threads, max_IJK + 1);
  cc_vec.allocate(n_threads, max_IJK + 1);
  eh_vec.allocate(n_threads, max_IJK + 1);

  E_s.initialise_fftw_plan(n_threads, eh_vec);
  H_s.initialise_fftw_plan(n_threads, eh_vec);

  N_e_x = I_tot - 1 + 1;
  N_e_y = J_tot - 1 + 1;
  N_e_z = K_tot - 1 + 1;
  N_h_x = I_tot + 1;
  N_h_y = J_tot + 1;
  N_h_z = K_tot + 1;

  dk_e_x = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_e_x));
  dk_e_y = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_e_y));
  dk_e_z = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_e_z));

  dk_h_x = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_h_x));
  dk_h_y = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_h_y));
  dk_h_z = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N_h_z));

  init_diff_shift_op(-0.5, dk_e_x, N_e_x);
  init_diff_shift_op(-0.5, dk_e_y, N_e_y);
  init_diff_shift_op(-0.5, dk_e_z, N_e_z);

  init_diff_shift_op(0.5, dk_h_x, N_h_x);
  init_diff_shift_op(0.5, dk_h_y, N_h_y);
  init_diff_shift_op(0.5, dk_h_z, N_h_z);
#endif

  params.set_Np(f_ex_vec);

  //fprintf(stderr,"Pre 01\n");
  //initialise E_norm and H_norm
  auto E_norm = (complex<double> *) malloc(f_ex_vec.size() * sizeof(complex<double>));
  auto H_norm = (complex<double> *) malloc(f_ex_vec.size() * sizeof(complex<double>));
  for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
    E_norm[ifx] = 0.;
    H_norm[ifx] = 0.;
  }

  //fprintf(stderr,"Pre 02\n");

  //fprintf(stderr,"Qos 00 (%d) (%d,%d,%d,%d):\n",J_tot,cuboid[0], cuboid[1],cuboid[4], cuboid[5]);
  /*set up surface mesh if required*/

  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    if (J_tot == 0)
      conciseCreateBoundary(cuboid[0], cuboid[1], cuboid[4], cuboid[5], &mx_surface_vertices,
                            &mx_surface_facets);
    else
      conciseTriangulateCuboidSkip(cuboid[0], cuboid[1], cuboid[2], cuboid[3], cuboid[4], cuboid[5],
                                   params.phasorinc, &mx_surface_vertices,
                                   &mx_surface_facets);
    //fprintf(stderr,"Qos 00a:\n");
    //we don't need the facets so destroy the matrix now to save memory
    mxDestroyArray(mx_surface_facets);
    dimptr_out = mxGetDimensions(mx_surface_vertices);
    n_surface_vertices = dimptr_out[0];
    //cast the vertex array as a 2-d integer array
    surface_vertices = cast_matlab_2D_array((int *) mxGetPr((mxArray *) mx_surface_vertices),
                                            dimptr_out[0], dimptr_out[1]);
    //create space for the complex amplitudes E and H around the surface. These will be in a large complex
    //array with each line being of the form Re(Ex) Im(Ex) Re(Ey) ... Im(Hz). Each line corresponds to the
    //the vertex with the same line as in surface_vertices
    ndims = 3;

    dims[0] = n_surface_vertices;
    dims[1] = 6;//one for each component of field
    dims[2] = f_ex_vec.size();

    mx_surface_amplitudes =
            mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    surface_EHr = cast_matlab_3D_array(mxGetPr((mxArray *) mx_surface_amplitudes), dims[0], dims[1],
                                       dims[2]);
    surface_EHi = cast_matlab_3D_array(mxGetPi((mxArray *) mx_surface_amplitudes), dims[0], dims[1],
                                       dims[2]);
    //now need to add a command to update the complex amplitudes
  }

  //fprintf(stderr,"Pre 03\n");
  /*Now set up the phasor array, we will have 3 complex output arrays for Ex, Ey and Ez.
    Phasors are extracted over the range Dxl + 3 - 1 to I_tot - Dxu - 1 to avoid pml cells
    see page III.80 for explanation of the following. This has been extended so that interpolation
    is done at the end of the FDTD run and also to handle the case of when there is no PML in place
    more appropriatley*/
  E.il = H.il = (params.pml.Dxl) ? params.pml.Dxl + 2 : 0;
  E.iu = H.iu = (params.pml.Dxu) ? I_tot - params.pml.Dxu - 1 : I_tot;
  E.jl = H.jl = (params.pml.Dyl) ? params.pml.Dyl + 2 : 0;
  E.ju = H.ju = (params.pml.Dyu) ? J_tot - params.pml.Dyu - 1 : J_tot;
  E.kl = H.kl = (params.pml.Dzl) ? params.pml.Dzl + 2 : 0;
  E.ku = H.ku = (params.pml.Dzu) ? K_tot - params.pml.Dzu - 1 : K_tot;

  E.I_tot = H.I_tot = E.iu - E.il + 1;
  E.J_tot = H.J_tot = E.ju - E.jl + 1;
  E.K_tot = H.K_tot = E.ku - E.kl + 1;

  //fprintf(stderr,"Pre 04\n");
  //fprintf(stderr,"pind_ju: %d, pind_jl: %d, J_tot: %d\n",pind_ju,pind_jl,J_tot);

  //fprintf(stderr,"Qos 01:\n");

  /*  dims[0] = I_tot - params.pml.Dxu - params.pml.Dxl - 3 + 1;
      dims[1] = J_tot - params.pml.Dyu - params.pml.Dyl - 3 + 1;
      dims[2] = K_tot - params.pml.Dzu - params.pml.Dzl - 3 + 1;*/
  auto output_grid_labels = GridLabels();

  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    ndims = 3;

    dims[0] = E.I_tot;
    dims[1] = E.J_tot;
    dims[2] = E.K_tot;

    fprintf(stderr, "dims:(%d,%d,%d)\n", dims[0], dims[1], dims[2]);

    plhs[0] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
    plhs[1] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
    plhs[2] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez

    plhs[3] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hx
    plhs[4] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hy
    plhs[5] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hz

    E.real.x = cast_matlab_3D_array(mxGetPr((mxArray *) plhs[0]), E.I_tot, E.J_tot, E.K_tot);
    E.imag.x = cast_matlab_3D_array(mxGetPi((mxArray *) plhs[0]), E.I_tot, E.J_tot, E.K_tot);

    E.real.y = cast_matlab_3D_array(mxGetPr((mxArray *) plhs[1]), E.I_tot, E.J_tot, E.K_tot);
    E.imag.y = cast_matlab_3D_array(mxGetPi((mxArray *) plhs[1]), E.I_tot, E.J_tot, E.K_tot);

    E.real.z = cast_matlab_3D_array(mxGetPr((mxArray *) plhs[2]), E.I_tot, E.J_tot, E.K_tot);
    E.imag.z = cast_matlab_3D_array(mxGetPi((mxArray *) plhs[2]), E.I_tot, E.J_tot, E.K_tot);

    H.real.x = cast_matlab_3D_array(mxGetPr((mxArray *) plhs[3]), H.I_tot, H.J_tot, H.K_tot);
    H.imag.x = cast_matlab_3D_array(mxGetPi((mxArray *) plhs[3]), H.I_tot, H.J_tot, H.K_tot);

    H.real.y = cast_matlab_3D_array(mxGetPr((mxArray *) plhs[4]), H.I_tot, H.J_tot, H.K_tot);
    H.imag.y = cast_matlab_3D_array(mxGetPi((mxArray *) plhs[4]), H.I_tot, H.J_tot, H.K_tot);

    H.real.z = cast_matlab_3D_array(mxGetPr((mxArray *) plhs[5]), H.I_tot, H.J_tot, H.K_tot);
    H.imag.z = cast_matlab_3D_array(mxGetPi((mxArray *) plhs[5]), H.I_tot, H.J_tot, H.K_tot);
    //fprintf(stderr,"Pre 05\n");
    //fprintf(stderr,"Qos 02:\n");
    //these will ultimately be copies of the phasors used to test convergence
    if (params.source_mode == SourceMode::steadystate) {
      dummy_array[0] =
              mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
      dummy_array[1] =
              mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
      dummy_array[2] =
              mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez
    }
    //fprintf(stderr,"Pre 06\n");
    //fprintf(stderr,"Qos 03:\n");
    if (params.source_mode == SourceMode::steadystate) {
      E_copy.real.x =
              cast_matlab_3D_array(mxGetPr((mxArray *) dummy_array[0]), dims[0], dims[1], dims[2]);
      E_copy.imag.x =
              cast_matlab_3D_array(mxGetPi((mxArray *) dummy_array[0]), dims[0], dims[1], dims[2]);

      E_copy.real.y =
              cast_matlab_3D_array(mxGetPr((mxArray *) dummy_array[1]), dims[0], dims[1], dims[2]);
      E_copy.imag.y =
              cast_matlab_3D_array(mxGetPi((mxArray *) dummy_array[1]), dims[0], dims[1], dims[2]);

      E_copy.real.z =
              cast_matlab_3D_array(mxGetPr((mxArray *) dummy_array[2]), dims[0], dims[1], dims[2]);
      E_copy.imag.z =
              cast_matlab_3D_array(mxGetPi((mxArray *) dummy_array[2]), dims[0], dims[1], dims[2]);
    }
    //fprintf(stderr,"Pre 07\n");
    //this will be a copy of the phasors which are extracted from the previous cycle

    //fprintf(stderr,"Qos 04:\n");

    //now construct the grid labels
    label_dims[0] = 1;
    label_dims[1] = dims[0];
    plhs[10] = mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//x
    output_grid_labels.x = mxGetPr((mxArray *) plhs[10]);

    label_dims[0] = 1;
    label_dims[1] = dims[1];
    //fprintf(stderr,"plhs[11]: %d,%d\n",label_dims[0],label_dims[1] );
    plhs[11] = mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//y
    output_grid_labels.y = mxGetPr((mxArray *) plhs[11]);

    label_dims[0] = 1;
    label_dims[1] = dims[2];
    plhs[12] = mxCreateNumericArray(2, (const mwSize *) label_dims, mxDOUBLE_CLASS, mxREAL);//y
    output_grid_labels.z = mxGetPr((mxArray *) plhs[12]);
    //fprintf(stderr,"Pre 08\n");
  } else {
    //initialise to empty matrices
    //fprintf(stderr,"Pre 09\n");
    ndims = 2;

    dims[0] = 0;
    dims[1] = 0;

    plhs[0] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ex
    plhs[1] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ey
    plhs[2] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Ez

    plhs[3] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hx
    plhs[4] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hy
    plhs[5] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);//Hz

    plhs[10] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                    mxCOMPLEX);//x_grid_labels_out
    plhs[11] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                    mxCOMPLEX);//y_grid_labels_out
    plhs[12] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS,
                                    mxCOMPLEX);//z_grid_labels_out
  }
  //fprintf(stderr,"Pre 10\n");
  //fprintf(stderr,"Qos 05:\n");
  //plhs[13] -> plhs[15] are the interpolated electric field values
  //plhs[16] -> plhs[18] are the interpolated magnetic field values

  //initialise arrays
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    E.zero();
    H.zero();
  }
  //fprintf(stderr,"Pre 11\n");
  if (params.exdetintegral && params.run_mode == RunMode::complete) {
    ndims = 2;
    dims[0] = 1;
    dims[1] = 1;
    const char *fieldnames[] = {"Idx", "Idy"};
    plhs[26] = mxCreateStructArray(ndims, (const mwSize *) dims, 2, fieldnames);

    ndims = 2;
    dims[0] = D_tilde.num_det_modes();
    dims[1] = f_ex_vec.size();

    mx_Idx = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    Idx_re = cast_matlab_2D_array(mxGetPr(mx_Idx), dims[0], dims[1]);
    Idx_im = cast_matlab_2D_array(mxGetPi(mx_Idx), dims[0], dims[1]);

    mx_Idy = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    Idy_re = cast_matlab_2D_array(mxGetPr(mx_Idy), dims[0], dims[1]);
    Idy_im = cast_matlab_2D_array(mxGetPi(mx_Idy), dims[0], dims[1]);

    Idx = (complex<double> **) malloc(sizeof(complex<double> *) * f_ex_vec.size());
    Idy = (complex<double> **) malloc(sizeof(complex<double> *) * f_ex_vec.size());

    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      Idx[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
      Idy[ifx] = (complex<double> *) malloc(sizeof(complex<double>) * dims[0]);
      for (int im = 0; im < dims[0]; im++) {
        Idx[ifx][im] = 0.;
        Idy[ifx][im] = 0.;
      }
    }

    for (int im = 0; im < dims[0]; im++) {
      for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
        Idx_re[ifx][im] = 0.;
        Idx_im[ifx][im] = 0.;
        Idy_re[ifx][im] = 0.;
        Idy_im[ifx][im] = 0.;
      }
    }

    mxSetField(plhs[26], 0, "Idx", mx_Idx);
    mxSetField(plhs[26], 0, "Idy", mx_Idy);
  } else {
    ndims = 2;
    dims[0] = 0;
    dims[1] = 0;

    plhs[26] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  }
  //fprintf(stderr,"Pre 12\n");
  if (params.run_mode == RunMode::complete && params.source_mode == SourceMode::steadystate && params.exphasorsvolume) {
    E_copy.zero();
  }
  //fprintf(stderr,"Pre 13\n");
  /*This is just for efficiency */
  K = K_tot - params.pml.Dxl - params.pml.Dxu;

  /*Now set up the phasor arrays for storing the fdtd version of the input fields,
    these will be used in a boot strapping procedure. Calculated over a complete
    xy-plane. */

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
  //work out if we have conductive background
  bool is_cond = is_conductive(rho_cond, I_tot, J_tot, K_tot);
  //work out if we have a dispersive background
  if (params.is_disp_ml) params.is_disp_ml = is_dispersive_ml(ml, K_tot);
  //  fprintf(stderr,"is_disp:%d, is_cond%d, params.is_disp_ml: %d\n",is_disp,is_cond,params.is_disp_ml);
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
  if (is_cond) { J_c.allocate_and_zero(); }
  /*end dispersive*/

  plhs[27] = fieldsample.mx;

  if (campssample.n_vertices() > 0) {
    ndims = 3;
    dims[0] = campssample.n_vertices();
    dims[1] = campssample.components.size();
    dims[2] = f_ex_vec.size();
    mx_camplitudes = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
    camplitudesR = cast_matlab_3D_array(mxGetPr(mx_camplitudes), dims[0], dims[1], dims[2]);
    camplitudesI = cast_matlab_3D_array(mxGetPi(mx_camplitudes), dims[0], dims[1], dims[2]);

  } else {
    ndims = 3;
    dims[0] = 0;
    dims[1] = 0;
    dims[2] = 0;
    mx_camplitudes = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  }
  plhs[28] = mx_camplitudes;
  //fprintf(stderr,"Pre 15\n");
  /*end of setup the output array for the sampled field*/


  /*set up the parameters for the phasor convergence procedure*/
  /*First we set dt so that an integer number of time periods fits within a sinusoidal period
   */
  double Nsteps_tmp = 0.0;
  double dt_old;
  if (params.source_mode == SourceMode::steadystate) {
    dt_old = params.dt;
    Nsteps_tmp = ceil(2. * dcpi / params.omega_an / params.dt * 3);
    params.dt = 2. * dcpi / params.omega_an * 3 / Nsteps_tmp;
  }

  //fprintf(stderr,"Pre 16\n");
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete)
    fprintf(stderr, "Changing dt from %.10e to %.10e\n", dt_old, params.dt);
  Nsteps = (int) lround(Nsteps_tmp);
  //fprintf(stderr,"Pre 17\n");
  //Nsteps = (int)(floor(3*2.*dcpi/(params.omega_an*params.dt)) + 1.);//the number of time steps in a sinusoidal period
  dft_counter = 0;
  //fprintf(stderr,"Pre 18\n");
  /*params.Nt should be an integer number of Nsteps in the case of steady-state operation*/
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete)
    if (params.Nt / Nsteps * Nsteps != params.Nt) {
      fprintf(stderr, "Changing the value of Nt from %d to", params.Nt);
      params.Nt = params.Nt / Nsteps * Nsteps;
      fprintf(stderr, " %d for correct phasor extraction\n", params.Nt);
    }
  //fprintf(stderr,"Pre 19\n");

  if ((params.run_mode == RunMode::complete) && (params.source_mode == SourceMode::steadystate)) printf("Nsteps: %d \n", Nsteps);

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
  //    fprintf(stderr,"params.start_tind: %d\n",params.start_tind);
  //  fprintf(stdout,"dz: %e, c: %e, dz/c: %e\n",dz,light_v,dz/light_v);
  //Begin of main iteration loop
  auto main_loop_timer = Timer();

  if (TIME_MAIN_LOOP) { main_loop_timer.start(); }

  for (unsigned int tind = params.start_tind; tind < params.Nt; tind++) {
    //fprintf(stderr,"Pos 00:\n");
    time_E = ((double) (tind + 1)) * params.dt;
    time_H = time_E - params.dt / 2.;
    //Extract phasors
    auto timer = Timer();
    if ((dft_counter == Nsteps) && (params.run_mode == RunMode::complete) && (params.source_mode == SourceMode::steadystate) &&
        params.exphasorsvolume) {//params.run_mode=complete,sourcemode=steadystate
      dft_counter = 0;
      double tol = checkPhasorConvergence(E, E_copy);

      //      mexPrintf("tol: %.5e \n",tol);
      fprintf(stderr, "tol: %.5e \n", tol);

      if (tol < TOL) break;//required accuracy obtained

      copyPhasors(E, E_copy, (int) mxGetNumberOfElements((mxArray *) plhs[0]));

      //clean the phasors
      E.zero();
      H.zero();

      if (params.exphasorssurface) {
        initialiseDouble3DArray(surface_EHr, n_surface_vertices, 6, f_ex_vec.size());
        initialiseDouble3DArray(surface_EHi, n_surface_vertices, 6, f_ex_vec.size());
      }
      //cleanphasors
    }
    //fprintf(stderr,"Pos 01:\n");

    if ((params.source_mode == SourceMode::steadystate) && (params.run_mode == RunMode::complete) && params.exphasorsvolume) {

      E.set_phasors(E_s, dft_counter - 1, params.omega_an, params.dt, Nsteps);
      H.set_phasors(H_s, dft_counter, params.omega_an, params.dt, Nsteps);

      if (params.exphasorssurface) {
        if (params.intphasorssurface) {
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++)
            extractPhasorsSurface(surface_EHr[ifx], surface_EHi[ifx], E_s, H_s, surface_vertices,
                                  n_surface_vertices, dft_counter, f_ex_vec[ifx] * 2 * dcpi,
                                  Nsteps, J_tot, params);
          dft_counter++;
        } else {
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++)
            extractPhasorsSurfaceNoInterpolation(surface_EHr[ifx], surface_EHi[ifx], E_s, H_s,
                                                 surface_vertices, n_surface_vertices, dft_counter,
                                                 f_ex_vec[ifx] * 2 * dcpi, Nsteps,
                                                 J_tot, params);
          dft_counter++;
        }
      }

    } else if ((params.source_mode == SourceMode::pulsed) && (params.run_mode == RunMode::complete) && params.exphasorsvolume) {
      if (TIME_EXEC) { timer.click(); }

      if ((tind - params.start_tind) % params.Np == 0) {
        E.set_phasors(E_s, tind - 1, params.omega_an, params.dt, params.Npe);
        H.set_phasors(H_s, tind, params.omega_an, params.dt, params.Npe);
      }
      if (TIME_EXEC) { timer.click(); }
      //fprintf(stderr,"Pos 01b:\n");
    }
    /*extract fieldsample*/
    if (fieldsample.all_vectors_are_non_empty()) {
      //if( (tind-params.start_tind) % params.Np == 0){
      double Ex_temp = 0., Ey_temp = 0., Ez_temp = 0.;

#pragma omp parallel default(shared) private(Ex_temp, Ey_temp, Ez_temp)
        {
#pragma omp for
          for (int kt = 0; kt < fieldsample.k.size(); kt++)
            for (int jt = 0; jt < fieldsample.j.size(); jt++)
              for (int it = 0; it < fieldsample.i.size(); it++) {
                ////fprintf(stderr,"Pos fs 1\n");
                interpolateTimeDomainFieldCentralEBandLimited(
                        E_s.xy, E_s.xz, E_s.yx, E_s.yz, E_s.zx, E_s.zy,
                        fieldsample.i[it] + params.pml.Dxl - 1,
                        fieldsample.j[jt] + params.pml.Dyl - 1,
                        fieldsample.k[kt] + params.pml.Dzl - 1, &Ex_temp, &Ey_temp, &Ez_temp);
                //fprintf(stderr,"Pos fs 2\n");
                for (int nt = 0; nt < fieldsample.n.size(); nt++)
                  fieldsample[nt][kt][jt][it] =
                          fieldsample[nt][kt][jt][it] +
                          pow(Ex_temp * Ex_temp + Ey_temp * Ey_temp + Ez_temp * Ez_temp,
                              fieldsample.n[nt] / 2.) /
                                  params.Nt;
                //fprintf(stderr,"%d %d %d %d -> %d %d %d (%d) %d [%d %d]\n",nt,kt,jt,it,(int)fieldsample_vecs.n[nt], (int)fieldsample_vecs.i[it] + params.pml.Dxl - 1, (int)fieldsample_vecs.j[jt] + params.pml.Dyl - 1, params.pml.Dyl,(int)fieldsample_vecs.k[kt] + params.pml.Dzl - 1 , Nsteps, (int)fieldsample_vecs.n[nt] - 2);
              }
        }
    }

    /*end extract fieldsample*/

    //fprintf(stderr,"Pos 02:\n");
    if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete && params.exphasorssurface) {
      if ((tind - params.start_tind) % params.Np == 0) {
        if (params.intphasorssurface)
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++)
            extractPhasorsSurface(surface_EHr[ifx], surface_EHi[ifx], E_s, H_s, surface_vertices,
                                  n_surface_vertices, tind, f_ex_vec[ifx] * 2 * dcpi, params.Npe, J_tot, params);
        else
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++)
            extractPhasorsSurfaceNoInterpolation(
                    surface_EHr[ifx], surface_EHi[ifx], E_s, H_s, surface_vertices,
                    n_surface_vertices, tind, f_ex_vec[ifx] * 2 * dcpi, params.Npe, J_tot, params);
      }
    }

    if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete && (campssample.n_vertices() > 0)) {
      //     fprintf(stderr,"loc 01 (%d,%d,%d)\n",tind,params.start_tind,params.Np);
      if ((tind - params.start_tind) % params.Np == 0) {
        //	fprintf(stderr,"loc 02\n");
        if (campssample.n_vertices() > 0) {
          //fprintf(stderr,"loc 03\n");
          //	  fprintf(stderr,"EPV 01\n");
          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            extractPhasorsVertices(camplitudesR[ifx], camplitudesI[ifx], E_s, H_s, campssample,
                                   tind, f_ex_vec[ifx] * 2 * dcpi, params.dt, params.Npe, params.dimension,
                                   J_tot, params.interp_method);
          }
        }
      }
    }
    
    //fprintf(stderr,"Pos 02a:\n");
    if (params.source_mode == SourceMode::pulsed && params.run_mode == RunMode::complete && params.exdetintegral) {
      if ((tind - params.start_tind) % params.Np == 0) {
        //First need to sum up the Ex and Ey values on a plane ready for FFT, remember that Ex_t and Ey_t are in row-major format whilst Exy etc. are in column major format
        for (j = params.pml.Dyl; j < (J_tot - params.pml.Dyu); j++)
          for (i = params.pml.Dxl; i < (I_tot - params.pml.Dxu); i++) {
            int n = j - params.pml.Dyl + (i - params.pml.Dxl) * (J_tot - params.pml.Dyu - params.pml.Dyl);
            Ex_t.v[n][0] = E_s.xy[params.k_det_obs][j][i] + E_s.xz[params.k_det_obs][j][i];
            Ex_t.v[n][1] = 0.;
            Ey_t.v[n][0] = E_s.yx[params.k_det_obs][j][i] + E_s.yz[params.k_det_obs][j][i];
            Ey_t.v[n][1] = 0.;
          }
        //fprintf(stderr,"Pos 02a [1] (%d,%d,%d,%d):\n",params.pml.Dyl,J_tot-params.pml.Dyu,params.pml.Dxl,I_tot-params.pml.Dxu);
        fftw_execute(Ex_t.plan);
        fftw_execute(Ey_t.plan);
        //fprintf(stderr,"Pos 02a [2]:\n");
        //Iterate over each mode
        for (int im = 0; im < D_tilde.num_det_modes(); im++) {
          //Now go back to column-major
          for (j = 0; j < (J_tot - params.pml.Dyu - params.pml.Dyl); j++)
            for (i = 0; i < (I_tot - params.pml.Dxu - params.pml.Dxl); i++) {
              int n = j + i * (J_tot - params.pml.Dyu - params.pml.Dyl);
              Ex_t.cm[j][i] = Ex_t.v[n][0] + I * Ex_t.v[n][1];
              Ey_t.cm[j][i] = Ey_t.v[n][0] + I * Ey_t.v[n][1];
            }
          //fprintf(stderr,"Pos 02a [3]:\n");
          //Now multiply the pupil, mostly the pupil is non-zero in only a elements
          for (j = 0; j < (J_tot - params.pml.Dyu - params.pml.Dyl); j++)
            for (i = 0; i < (I_tot - params.pml.Dxu - params.pml.Dxl); i++) {
              Ex_t.cm[j][i] = Ex_t.cm[j][i] * pupil[j][i] * D_tilde.x[j][i][im];
              Ey_t.cm[j][i] = Ey_t.cm[j][i] * pupil[j][i] * D_tilde.y[j][i][im];
            }
            //fprintf(stderr,"Pos 02a [4]:\n");
            //now iterate over each frequency to extract phasors at
#pragma omp parallel default(shared) private(lambda_an_t, Idxt, Idyt, i, j, kprop, phaseTermE,     \
                                             cphaseTermE)
          {
#pragma omp for
            for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
              //wavelength in air
              lambda_an_t = light_v / f_ex_vec[ifx];
              //fprintf(stdout,"lambda_an_t = %e, light_v = %e, z_obs = %e\n",lambda_an_t,light_v,z_obs);
              Idxt = 0.;
              Idyt = 0.;

              //now loop over all angular frequencies
              for (j = 0; j < (J_tot - params.pml.Dyu - params.pml.Dyl); j++)
                for (i = 0; i < (I_tot - params.pml.Dxu - params.pml.Dxl); i++) {
                  if ((lambda_an_t * f_vec.x[i] * lambda_an_t * f_vec.x[i] +
                       lambda_an_t * f_vec.y[j] * lambda_an_t * f_vec.y[j]) < 1) {

                    if (!params.air_interface_present) {
                      /*This had to be fixed since we must take into account the refractive index of the medium.

           */
                      kprop = exp(I * params.z_obs * 2. * dcpi / lambda_an_t * refind *
                                  sqrt(1. - pow(lambda_an_t * f_vec.x[i] / refind, 2.) -
                                       pow(lambda_an_t * f_vec.y[j] / refind, 2.)));
                      //fprintf(stdout,"%d %d %e %e %e %e %e %e %e\n",i,j,f_vec.x[i],f_vec.y[j],real(kprop),imag(kprop),z_obs,dcpi,lambda_an_t);
                    } else {
                      kprop = exp(I * (-params.air_interface + params.z_obs) * 2. * dcpi / lambda_an_t * refind *
                                  sqrt(1. - pow(lambda_an_t * f_vec.x[i] / refind, 2.) -
                                       pow(lambda_an_t * f_vec.y[j] / refind, 2.))) *
                              exp(I * params.air_interface * 2. * dcpi / lambda_an_t *
                                  sqrt(1. - pow(lambda_an_t * f_vec.x[i], 2.) -
                                       pow(lambda_an_t * f_vec.y[j], 2.)));
                    }
                  } else
                    kprop = 0.;

                  Idxt += Ex_t.cm[j][i] * kprop;
                  Idyt += Ey_t.cm[j][i] * kprop;
                }
              phaseTermE = fmod(f_ex_vec[ifx] * 2. * dcpi * ((double) tind) * params.dt, 2 * dcpi);
              cphaseTermE = exp(phaseTermE * I) * 1. / ((double) params.Npe);

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
                            I_tot, J_tot, K0.index + 1, tind, params.omega_an, params.dt,
                            params.Nt);//extract the phasors just above the line
      }
    //fprintf(stderr,"Pos 02c:\n");

    //Update equations for the E field

    /*There are two options for determing the update coefficients for the FDTD cell:

      1) If cell (i,j,k) is either free space or PML:

      materials[k][j][i] will be set to 0. In this case the update parameter used will
      be given by C.a.y[j], C.b.y[j] etc depending on which update equation is being implemented.

      2) if cell (i,j,k) is composed of a scattering type material then materials[k][j][i] will be
      non-zero and will be an index into Cmaterial.a.y and Cmaterial.b.y etc depending on which
      update equation is being implemented.

    */

    int array_ind = 0;
    //fprintf(stderr,"I_tot=%d, J_tot=%d, K_tot=%d\n",I_tot,J_tot,K_tot);
    if (TIME_EXEC) { timer.click(); }
    //fprintf(stderr,"Dimension = %d\n",params.dimension);
    /*
      for(k=0;k<(K_tot+1);k++)
      fprintf(stdout,"%e ",Exy[k][13][13]+Exz[k][13][13]);
      fprintf(stdout,"\n");
    */
#pragma omp parallel default(shared) private(i, j, k, n, rho, k_loc, array_ind, Ca, Cb, Cc, alpha_l,\
                                             beta_l, gamma_l, kappa_l, sigma_l, Enp1,              \
                                             Jnp1)//,ca_vec,cb_vec,cc_vec,eh_vec)
    {
      n = omp_get_thread_num();
      Enp1 = 0.0;
      array_ind = 0;

      if (params.dimension == THREE || params.dimension == TE) {
#ifdef FDFLAG// Use central difference derivatives
             //FDTD, Exy
#pragma omp for
        for (k = 0; k < (K_tot + 1); k++)
          for (j = 1; j < J_tot; j++)
            for (i = 0; i < I_tot; i++) {
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
                if (is_cond) rho = rho_cond.y[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.y[array_ind];
                kappa_l = ml.kappa.y[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][j][i + 1]) {
                    alpha_l += alpha[materials[k][j][i + 1] - 1];
                    beta_l += beta[materials[k][j][i + 1] - 1];
                    gamma_l += gamma[materials[k][j][i + 1] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }


              Enp1 = Ca * Exy[k][j][i] +
                     Cb * (Hzy[k][j][i] + Hzx[k][j][i] - Hzy[k][j - 1][i] - Hzx[k][j - 1][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.xy[k][j][i] -
                        1. / 2. * Cb * dy *
                                ((1 + alpha_l) * Jxy[k][j][i] + beta_l * J_nm1.xy[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dy * J_c.xy[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * Jxy[k][j][i] + beta_l * J_nm1.xy[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xy[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * Exy[k][j][i];

                E_nm1.xy[k][j][i] = Exy[k][j][i];
                J_nm1.xy[k][j][i] = Jxy[k][j][i];
                Jxy[k][j][i] = Jnp1;

                //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,Jxy[k][j][i]);
              }

              if (is_cond && rho) { J_c.xy[k][j][i] -= rho * (Enp1 + Exy[k][j][i]); }

              Exy[k][j][i] = Enp1;
            }
            //FDTD, Exy
#else//PSTD, Exy
//fprintf(stderr,"Pos 02d:\n");
#pragma omp for
        for (k = 0; k < (K_tot + 1); k++)
          for (i = 0; i < I_tot; i++) {
            for (j = 1; j < J_tot; j++) {
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
                if (is_cond) rho = rho_cond.y[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.y[array_ind];
                kappa_l = ml.kappa.y[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][j][i + 1]) {
                    alpha_l += alpha[materials[k][j][i + 1] - 1];
                    beta_l += beta[materials[k][j][i + 1] - 1];
                    gamma_l += gamma[materials[k][j][i + 1] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }


              Enp1 = 0.0;
              //Enp1 = Ca*Exy[k][j][i]+Cb*(Hzy[k][j][i] + Hzx[k][j][i] - Hzy[k][j-1][i] - Hzx[k][j-1][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.xy[k][j][i] -
                        1. / 2. * Cb * dy *
                                ((1 + alpha_l) * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dy * J_c.xy[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.xy[k][j][i] + beta_l * J_nm1.xy[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xy[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * E_s.xy[k][j][i];

                E_nm1.xy[k][j][i] = E_s.xy[k][j][i];
                J_nm1.xy[k][j][i] = J_s.xy[k][j][i];
                J_s.xy[k][j][i] = Jnp1;

                //	    fprintf(stderr,"(%d,%d,%d): %e\n",i,j,k,Jxy[k][j][i]);
              }

              if (is_cond && rho) { J_c.xy[k][j][i] -= rho * (Enp1 + E_s.xy[k][j][i]); }

              eh_vec[n][j][0] = H_s.zy[k][j][i] + H_s.zx[k][j][i];
              eh_vec[n][j][1] = 0.;
              ca_vec[n][j - 1] = Ca;
              cb_vec[n][j - 1] = Cb;
            }
            if (J_tot > 1) {
              j = 0;
              eh_vec[n][j][0] = H_s.zy[k][j][i] + H_s.zx[k][j][i];
              eh_vec[n][j][1] = 0.;
              first_derivative(eh_vec[n], eh_vec[n], dk_e_y,
                               N_e_y, E_s.xy.plan_f[n], E_s.xy.plan_b[n]);


              //fprintf(stdout,"(%d,%d) %d (of %d)\n",i,k,n,omp_get_num_threads());

              for (j = 1; j < J_tot; j++) {
                E_s.xy[k][j][i] = ca_vec[n][j - 1] * E_s.xy[k][j][i] +
                                  cb_vec[n][j - 1] *
                                          eh_vec[n][j][0] / ((double) N_e_y);
              }
            }
          }
          //PSTD, Exy
#endif

            /*
    if(is_disp){
    i=36;
    j=36;
    k=36;

    fprintf(stdout,"%e %e",Jxy[k][j][i],Exy[k][j][i]);
    }
  */

            //fprintf(stderr,"Pos 04:\n");
            //Exz updates
#ifdef FDFLAG// Use central difference derivatives
#pragma omp for
        for (k = 1; k < K_tot; k++)
          for (j = 0; j < J_tot_p1_bound; j++)
            for (i = 0; i < I_tot; i++) {
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
                if (is_cond) rho = rho_cond.z[k_loc];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.z[k_loc];
                kappa_l = ml.kappa.z[k_loc];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][j][i + 1]) {
                    alpha_l += alpha[materials[k][j][i + 1] - 1];
                    beta_l += beta[materials[k][j][i + 1] - 1];
                    gamma_l += gamma[materials[k][j][i + 1] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }
              /*if( materials[k][j][i] || materials[k][j][i+1])
      fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_cond:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_cond,rho,is_disp,params.is_disp_ml);
      if(tind==0)
      fprintf(stdout,"%d %d %e %e\n",i,k,Ca, Cb);*/
              Enp1 = Ca * Exz[k][j][i] +
                     Cb * (Hyx[k - 1][j][i] + Hyz[k - 1][j][i] - Hyx[k][j][i] - Hyz[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.xz[k][j][i] -
                        1. / 2. * Cb * dz *
                                ((1 + alpha_l) * Jxz[k][j][i] + beta_l * J_nm1.xz[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dz * J_c.xz[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * Jxz[k][j][i] + beta_l * J_nm1.xz[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xz[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * Exz[k][j][i];
                E_nm1.xz[k][j][i] = Exz[k][j][i];
                J_nm1.xz[k][j][i] = Jxz[k][j][i];
                Jxz[k][j][i] = Jnp1;
              }

              if (is_cond && rho) { J_c.xz[k][j][i] -= rho * (Enp1 + Exz[k][j][i]); }

              Exz[k][j][i] = Enp1;
            }
            //FDTD, Exz
#else//PSTD, Exz
        //#pragma omp for
        for (j = 0; j < J_tot_p1_bound; j++)
#pragma omp for
          for (i = 0; i < I_tot; i++) {
            for (k = 1; k < K_tot; k++) {
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
                if (is_cond) rho = rho_cond.z[k_loc];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.z[k_loc];
                kappa_l = ml.kappa.z[k_loc];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][j][i + 1]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][j][i + 1]) {
                    alpha_l += alpha[materials[k][j][i + 1] - 1];
                    beta_l += beta[materials[k][j][i + 1] - 1];
                    gamma_l += gamma[materials[k][j][i + 1] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }
              /*if( materials[k][j][i] || materials[k][j][i+1])
      fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_cond:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_cond,rho,is_disp,params.is_disp_ml);
      if(tind==0)
      fprintf(stdout,"%d %d %e %e\n",i,k,Ca, Cb);*/
              //Enp1 = Ca*Exz[k][j][i]+Cb*(Hyx[k-1][j][i] + Hyz[k-1][j][i] - Hyx[k][j][i] - Hyz[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.xz[k][j][i] -
                        1. / 2. * Cb * dz *
                                ((1 + alpha_l) * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dz * J_c.xz[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.xz[k][j][i] + beta_l * J_nm1.xz[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.xz[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * E_s.xz[k][j][i];
                E_nm1.xz[k][j][i] = E_s.xz[k][j][i];
                J_nm1.xz[k][j][i] = J_s.xz[k][j][i];
                J_s.xz[k][j][i] = Jnp1;
              }

              if (is_cond && rho) { J_c.xz[k][j][i] -= rho * (Enp1 + E_s.xz[k][j][i]); }

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
            first_derivative(eh_vec[n], eh_vec[n], dk_e_z,
                             N_e_z, E_s.xz.plan_f[n], E_s.xz.plan_b[n]);
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
                                cb_vec[n][k - 1] *
                                        eh_vec[n][k][0] / ((double) N_e_z);
            }
          }
          //PSTD, Exz
#endif

            //fprintf(stderr,"Pos 05:\n");
            //Eyx updates
#ifdef FDFLAG// Use central difference derivatives
             //FDTD, Eyx
#pragma omp for
        for (k = 0; k < (K_tot + 1); k++)
          for (j = 0; j < J_tot_bound; j++)
            for (i = 1; i < I_tot; i++) {
              rho = 0.;
              k_loc = k;
              if (params.is_structure) {
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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
                if (is_cond) rho = rho_cond.x[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.x[array_ind];
                kappa_l = ml.kappa.x[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][min(J_tot, j + 1)][i]) {
                    alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                    beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                    gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }


              Enp1 = Ca * Eyx[k][j][i] +
                     Cb * (Hzx[k][j][i - 1] + Hzy[k][j][i - 1] - Hzx[k][j][i] - Hzy[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.yx[k][j][i] -
                        1. / 2. * Cb * dx *
                                ((1 + alpha_l) * Jyx[k][j][i] + beta_l * J_nm1.yx[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dx * J_c.yx[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * Jyx[k][j][i] + beta_l * J_nm1.yx[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yx[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * Eyx[k][j][i];
                E_nm1.yx[k][j][i] = Eyx[k][j][i];
                J_nm1.yx[k][j][i] = Jyx[k][j][i];
                Jyx[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.yx[k][j][i] -= rho * (Enp1 + Eyx[k][j][i]); }

              Eyx[k][j][i] = Enp1;
            }
            //FDTD, Eyx
#else//PSTD, Eyx
#pragma omp for
        for (k = 0; k < (K_tot + 1); k++)
          for (j = 0; j < J_tot_bound; j++) {
            for (i = 1; i < I_tot; i++) {
              rho = 0.;
              k_loc = k;
              if (params.is_structure) {
                if (k > params.pml.Dzl && k < (params.pml.Dzl + K)) {
                  if ((k - structure[i][1]) < (K + params.pml.Dzl) && (k - structure[i][1]) > params.pml.Dzl)
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
                if (is_cond) rho = rho_cond.x[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.x[array_ind];
                kappa_l = ml.kappa.x[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][min(J_tot, j + 1)][i]) {
                    alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                    beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                    gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }


              //Enp1 = Ca*Eyx[k][j][i]+Cb*(Hzx[k][j][i-1] + Hzy[k][j][i-1] - Hzx[k][j][i] - Hzy[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.yx[k][j][i] -
                        1. / 2. * Cb * dx *
                                ((1 + alpha_l) * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dx * J_c.yx[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.yx[k][j][i] + beta_l * J_nm1.yx[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yx[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * E_s.yx[k][j][i];
                E_nm1.yx[k][j][i] = E_s.yx[k][j][i];
                J_nm1.yx[k][j][i] = J_s.yx[k][j][i];
                J_s.yx[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.yx[k][j][i] -= rho * (Enp1 + E_s.yx[k][j][i]); }

              eh_vec[n][i][0] = H_s.zx[k][j][i] + H_s.zy[k][j][i];
              eh_vec[n][i][1] = 0.;
              ca_vec[n][i - 1] = Ca;
              cb_vec[n][i - 1] = Cb;
            }
            i = 0;
            eh_vec[n][i][0] = H_s.zx[k][j][i] + H_s.zy[k][j][i];
            eh_vec[n][i][1] = 0.;

            first_derivative(eh_vec[n], eh_vec[n], dk_e_x,
                             N_e_x, E_s.yx.plan_f[n], E_s.yx.plan_b[n]);

            for (i = 1; i < I_tot; i++) {
              E_s.yx[k][j][i] = ca_vec[n][i - 1] * E_s.yx[k][j][i] -
                                cb_vec[n][i - 1] *
                                        eh_vec[n][i][0] / ((double) N_e_x);
              //Eyx[k][j][i] = Enp1;
            }
          }
          //PSTD, Eyx
#endif

            //fprintf(stderr,"Pos 06:\n");
            //Eyz updates
#ifdef FDFLAG// Use central difference derivatives
//FDTD, Eyz
#pragma omp for
        for (k = 1; k < K_tot; k++)
          for (j = 0; j < J_tot_bound; j++)
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
                if (is_cond) rho = rho_cond.z[k_loc];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.z[k_loc];
                kappa_l = ml.kappa.z[k_loc];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][min(J_tot, j + 1)][i]) {
                    alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                    beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                    gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }

              //fprintf(stderr,"[%d %d %d]Ca: %e, Cb: %e, Cc: %e, alpha: %e, beta: %e, gamme: %e\n",i,j,k,Ca,Cb,Cc,alpha_l,beta_l,gamma_l);
              Enp1 = Ca * Eyz[k][j][i] +
                     Cb * (Hxy[k][j][i] + Hxz[k][j][i] - Hxy[k - 1][j][i] - Hxz[k - 1][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.yz[k][j][i] -
                        1. / 2. * Cb * dz *
                                ((1 + alpha_l) * Jyz[k][j][i] + beta_l * J_nm1.yz[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dz * J_c.yz[k][j][i];

              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * Jyz[k][j][i] + beta_l * J_nm1.yz[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yz[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * Eyz[k][j][i];
                E_nm1.yz[k][j][i] = Eyz[k][j][i];
                J_nm1.yz[k][j][i] = Jyz[k][j][i];
                Jyz[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.yz[k][j][i] -= rho * (Enp1 + Eyz[k][j][i]); }

              Eyz[k][j][i] = Enp1;
            }
//FDTD, Eyz
#else//PSTD, Eyz
#pragma omp for
        for (j = 0; j < J_tot_bound; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            for (k = 1; k < K_tot; k++) {
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
                if (is_cond) rho = rho_cond.z[k_loc];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.z[k_loc];
                kappa_l = ml.kappa.z[k_loc];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k][min(J_tot, j + 1)][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k][min(J_tot, j + 1)][i]) {
                    alpha_l += alpha[materials[k][min(J_tot, j + 1)][i] - 1];
                    beta_l += beta[materials[k][min(J_tot, j + 1)][i] - 1];
                    gamma_l += gamma[materials[k][min(J_tot, j + 1)][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }

              //fprintf(stderr,"[%d %d %d]Ca: %e, Cb: %e, Cc: %e, alpha: %e, beta: %e, gamme: %e\n",i,j,k,Ca,Cb,Cc,alpha_l,beta_l,gamma_l);
              //Enp1 = Ca*Eyz[k][j][i]+Cb*(Hxy[k][j][i] + Hxz[k][j][i] - Hxy[k-1][j][i] - Hxz[k-1][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.yz[k][j][i] -
                        1. / 2. * Cb * dz *
                                ((1 + alpha_l) * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dz * J_c.yz[k][j][i];

              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.yz[k][j][i] + beta_l * J_nm1.yz[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.yz[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * E_s.yz[k][j][i];
                E_nm1.yz[k][j][i] = E_s.yz[k][j][i];
                J_nm1.yz[k][j][i] = J_s.yz[k][j][i];
                J_s.yz[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.yz[k][j][i] -= rho * (Enp1 + E_s.yz[k][j][i]); }

              eh_vec[n][k][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
              eh_vec[n][k][1] = 0.;
              ca_vec[n][k - 1] = Ca;
              cb_vec[n][k - 1] = Cb;
            }
            k = 0;
            eh_vec[n][k][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
            eh_vec[n][k][1] = 0.;
            first_derivative(eh_vec[n], eh_vec[n], dk_e_z,
                             N_e_z, E_s.yz.plan_f[n], E_s.yz.plan_b[n]);


            for (k = 1; k < K_tot; k++) {
              E_s.yz[k][j][i] = ca_vec[n][k - 1] * E_s.yz[k][j][i] +
                                cb_vec[n][k - 1] *
                                        eh_vec[n][k][0] / ((double) N_e_z);
              //Eyz[k][j][i] = Enp1;
            }
          }
          //PSTD, Eyz
#endif
      }//if(params.dimension==THREE || params.dimension==TE)

      //fprintf(stderr,"Pos 07:\n");
      if (params.dimension == THREE || params.dimension == TE) {
#ifdef FDFLAG// Use central difference derivatives
#pragma omp for
        //Ezx updates
        for (k = 0; k < K_tot; k++)
          for (j = 0; j < J_tot_p1_bound; j++)
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
                if (is_cond) rho = rho_cond.x[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.x[array_ind];
                kappa_l = ml.kappa.x[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k + 1][j][i]) {
                    alpha_l += alpha[materials[k + 1][j][i] - 1];
                    beta_l += beta[materials[k + 1][j][i] - 1];
                    gamma_l += gamma[materials[k + 1][j][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }

                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }

              /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_cond:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_cond,rho,is_disp,params.is_disp_ml);*/
              Enp1 = Ca * Ezx[k][j][i] +
                     Cb * (Hyx[k][j][i] + Hyz[k][j][i] - Hyx[k][j][i - 1] - Hyz[k][j][i - 1]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zx[k][j][i] -
                        1. / 2. * Cb * dx *
                                ((1 + alpha_l) * Jzx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dx * J_c.zx[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * Jzx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * Ezx[k][j][i];
                E_nm1.zx[k][j][i] = Ezx[k][j][i];
                J_nm1.zx[k][j][i] = Jzx[k][j][i];
                Jzx[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + Ezx[k][j][i]); }

              Ezx[k][j][i] = Enp1;
            }
//FDTD, Ezx
#else//PSTD, Ezx
#pragma omp for
        for (k = 0; k < K_tot; k++)
          for (j = 0; j < J_tot_p1_bound; j++) {
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
                if (is_cond) rho = rho_cond.x[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.x[array_ind];
                kappa_l = ml.kappa.x[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k + 1][j][i]) {
                    alpha_l += alpha[materials[k + 1][j][i] - 1];
                    beta_l += beta[materials[k + 1][j][i] - 1];
                    gamma_l += gamma[materials[k + 1][j][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }

                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }

              /*if( materials[k][j][i] || materials[k][j][i+1])
        fprintf(stdout,"(%d,%d,%d), Ca= %e, Cb=%e, is_cond:%d, rho: %e, is_disp: %d, params.is_disp_ml: %d\n",i,j,k,Ca,Cb,is_cond,rho,is_disp,params.is_disp_ml);*/
              //Enp1 = Ca*Ezx[k][j][i]+Cb*(Hyx[k][j][i] + Hyz[k][j][i] - Hyx[k][j][i-1] - Hyz[k][j][i-1]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zx[k][j][i] -
                        1. / 2. * Cb * dx *
                                ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dx * J_c.zx[k][j][i];
              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * E_s.zx[k][j][i];
                E_nm1.zx[k][j][i] = E_s.zx[k][j][i];
                J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                J_s.zx[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + E_s.zx[k][j][i]); }

              eh_vec[n][i][0] = H_s.yx[k][j][i] + H_s.yz[k][j][i];
              eh_vec[n][i][1] = 0.;
              ca_vec[n][i - 1] = Ca;
              cb_vec[n][i - 1] = Cb;
            }
            i = 0;
            eh_vec[n][i][0] = H_s.yx[k][j][i] + H_s.yz[k][j][i];
            eh_vec[n][i][1] = 0.;

            first_derivative(eh_vec[n], eh_vec[n], dk_e_x,
                             N_e_x, E_s.zx.plan_f[n], E_s.zx.plan_b[n]);

            for (i = 1; i < I_tot; i++) {
              E_s.zx[k][j][i] = ca_vec[n][i - 1] * E_s.zx[k][j][i] +
                                cb_vec[n][i - 1] *
                                        eh_vec[n][i][0] / ((double) N_e_x);
              //Ezx[k][j][i] = Enp1;
            }
          }
          //PSTD, Ezx
#endif
      }//(params.dimension==THREE || params.dimension==TE)
      else {
#pragma omp for
        //Ezx updates
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
                if (is_cond) rho = rho_cond.x[i];
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
                sigma_l = ml.sigma.x[array_ind];
                kappa_l = ml.kappa.x[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];

                if (materials[k][j][i]) {
                  alpha_l = alpha[materials[k][j][i] - 1];
                  beta_l = beta[materials[k][j][i] - 1];
                  gamma_l = gamma[materials[k][j][i] - 1];

                } else {
                  alpha_l = ml.alpha[k_loc];
                  beta_l = ml.beta[k_loc];
                  gamma_l = ml.gamma[k_loc];
                }
              }

              Enp1 = Ca * E_s.zx[k][j][i] + Cb * (H_s.yx[k][j][i] + H_s.yz[k][j][i] -
                                                  H_s.yx[k][j][i - 1] - H_s.yz[k][j][i - 1]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zx[k][j][i] -
                        1. / 2. * Cb * dx *
                                ((1 + alpha_l) * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dx * J_c.zx[k][j][i];

              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zx[k][j][i] + beta_l * J_nm1.zx[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zx[k][j][i]);
                Jnp1 += sigma_l / eo * gamma_l * E_s.zx[k][j][i];
                E_nm1.zx[k][j][i] = E_s.zx[k][j][i];
                J_nm1.zx[k][j][i] = J_s.zx[k][j][i];
                J_s.zx[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.zx[k][j][i] -= rho * (Enp1 + E_s.zx[k][j][i]); }

              E_s.zx[k][j][i] = Enp1;
            }
      }
      //fprintf(stderr,"Pos 08:\n");
      if (params.dimension == THREE || params.dimension == TE) {
#ifdef FDFLAG// Use central difference derivatives
             //FDTD, Ezy
#pragma omp for
        //Ezy updates
        for (k = 0; k < K_tot; k++)
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
                if (is_cond) rho = rho_cond.y[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.y[array_ind];
                kappa_l = ml.kappa.y[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k + 1][j][i]) {
                    alpha_l += alpha[materials[k + 1][j][i] - 1];
                    beta_l += beta[materials[k + 1][j][i] - 1];
                    gamma_l += gamma[materials[k + 1][j][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }


              Enp1 = Ca * Ezy[k][j][i] +
                     Cb * (Hxy[k][j - 1][i] + Hxz[k][j - 1][i] - Hxy[k][j][i] - Hxz[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zy[k][j][i] -
                        1. / 2. * Cb * dy *
                                ((1 + alpha_l) * Jzy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dy * J_c.zy[k][j][i];

              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * Jzy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                Jnp1 += sigma_l / eo * gamma_l * Ezy[k][j][i];
                E_nm1.zy[k][j][i] = Ezy[k][j][i];
                J_nm1.zy[k][j][i] = Jzy[k][j][i];
                Jzy[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + Ezy[k][j][i]); }
              Ezy[k][j][i] = Enp1;
            }
//FDTD, Ezy
#else//PSTD, Ezy
#pragma omp for
        //Ezy updates
        for (k = 0; k < K_tot; k++)
          for (i = 0; i < (I_tot + 1); i++) {
            for (j = 1; j < J_tot; j++) {
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
                if (is_cond) rho = rho_cond.y[array_ind];
              }

              alpha_l = 0.;
              beta_l = 0.;
              gamma_l = 0.;
              kappa_l = 1.;
              sigma_l = 0.;

              if (is_disp || params.is_disp_ml) {
                sigma_l = ml.sigma.y[array_ind];
                kappa_l = ml.kappa.y[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];
                if (materials[k][j][i] || materials[k + 1][j][i]) {
                  if (materials[k][j][i]) {
                    alpha_l = alpha[materials[k][j][i] - 1];
                    beta_l = beta[materials[k][j][i] - 1];
                    gamma_l = gamma[materials[k][j][i] - 1];
                  } else {
                    alpha_l = ml.alpha[k_loc];
                    beta_l = ml.beta[k_loc];
                    gamma_l = ml.gamma[k_loc];
                  }

                  if (materials[k + 1][j][i]) {
                    alpha_l += alpha[materials[k + 1][j][i] - 1];
                    beta_l += beta[materials[k + 1][j][i] - 1];
                    gamma_l += gamma[materials[k + 1][j][i] - 1];
                  } else {
                    alpha_l += ml.alpha[k_loc];
                    beta_l += ml.beta[k_loc];
                    gamma_l += ml.gamma[k_loc];
                  }
                  alpha_l = alpha_l / 2.;
                  beta_l = beta_l / 2.;
                  gamma_l = gamma_l / 2.;
                }
              }


              //Enp1 = Ca*Ezy[k][j][i]+Cb*(Hxy[k][j-1][i] + Hxz[k][j-1][i] - Hxy[k][j][i] - Hxz[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zy[k][j][i] -
                        1. / 2. * Cb * dy *
                                ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dy * J_c.zy[k][j][i];

              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                Jnp1 += sigma_l / eo * gamma_l * E_s.zy[k][j][i];
                E_nm1.zy[k][j][i] = E_s.zy[k][j][i];
                J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                J_s.zy[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + E_s.zy[k][j][i]); }

              eh_vec[n][j][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
              eh_vec[n][j][1] = 0.;
              ca_vec[n][j - 1] = Ca;
              cb_vec[n][j - 1] = Cb;
            }
            if (J_tot > 1) {
              j = 0;
              eh_vec[n][j][0] = H_s.xy[k][j][i] + H_s.xz[k][j][i];
              eh_vec[n][j][1] = 0.;
              first_derivative(eh_vec[n], eh_vec[n], dk_e_y,
                               N_e_y, E_s.zy.plan_f[n], E_s.zy.plan_b[n]);
            }
            for (j = 1; j < J_tot; j++) {
              E_s.zy[k][j][i] = ca_vec[n][j - 1] * E_s.zy[k][j][i] -
                                cb_vec[n][j - 1] *
                                        eh_vec[n][j][0] / ((double) N_e_y);
              //Ezy[k][j][i] = Enp1;
            }
          }
//PSTD, Ezy
#endif
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
                if (is_cond) rho = rho_cond.y[array_ind];
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
                kappa_l = ml.kappa.y[array_ind];
                sigma_l = ml.sigma.y[array_ind];
                alpha_l = ml.alpha[k_loc];
                beta_l = ml.beta[k_loc];
                gamma_l = ml.gamma[k_loc];

                if (!materials[k][j][i]) {
                  alpha_l = 0.;
                  beta_l = 0.;
                  gamma_l = 0.;
                } else {
                  alpha_l = ml.alpha[k_loc];
                  beta_l = ml.beta[k_loc];
                  gamma_l = ml.gamma[k_loc];
                }
              }


              Enp1 = Ca * E_s.zy[k][j][i] + Cb * (H_s.xy[k][j - 1][i] + H_s.xz[k][j - 1][i] -
                                                  H_s.xy[k][j][i] - H_s.xz[k][j][i]);
              if ((is_disp || params.is_disp_ml) && gamma_l)
                Enp1 += Cc * E_nm1.zy[k][j][i] -
                        1. / 2. * Cb * dy *
                                ((1 + alpha_l) * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i]);
              if (is_cond && rho) Enp1 += Cb * dy * J_c.zy[k][j][i];

              if ((is_disp || params.is_disp_ml) && gamma_l) {
                Jnp1 = alpha_l * J_s.zy[k][j][i] + beta_l * J_nm1.zy[k][j][i] +
                       kappa_l * gamma_l / (2. * params.dt) * (Enp1 - E_nm1.zy[k][j][i]);

                Jnp1 += sigma_l / eo * gamma_l * E_s.zy[k][j][i];
                E_nm1.zy[k][j][i] = E_s.zy[k][j][i];
                J_nm1.zy[k][j][i] = J_s.zy[k][j][i];
                J_s.zy[k][j][i] = Jnp1;
              }
              if (is_cond && rho) { J_c.zy[k][j][i] -= rho * (Enp1 + E_s.zy[k][j][i]); }

              E_s.zy[k][j][i] = Enp1;
            }
      }
    }//end of parallel section
    //fprintf(stderr,"Pos 09:\n");
    if (TIME_EXEC) { timer.click(); }
    /********************/

    //update terms for self consistency across scattered/total interface - E updates##
    if (params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase = exp(-I * fmod(params.omega_an * time_H, 2. * dcpi));
      double commonAmplitude = linearRamp(time_H, 1. / (params.omega_an / (2 * dcpi)), ramp_width);
      for (k = (K0.index); k <= (K1.index); k++)
        for (j = (J0.index); j <= (J1.index); j++) {
          if (I0.apply) {//Perform across I0

            if (!params.is_multilayer) array_ind = I0.index;
            else
              array_ind = (I_tot + 1) * k + I0.index;

            if (k < (K1.index) || params.dimension == TM) {
              E_s.zx[k][j][I0.index] =
                      E_s.zx[k][j][I0.index] -
                      C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
              if (is_cond)
                J_c.zx[k][j][I0.index] +=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
              if (params.is_disp_ml)
                J_s.zx[k][j][I0.index] +=
                        ml.kappa.x[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][2] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][2]));
            }
            if (j < (J1.index)) {
              E_s.yx[k][j][I0.index] =
                      E_s.yx[k][j][I0.index] +
                      C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
              if (is_cond)
                J_c.yx[k][j][I0.index] -=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
              if (params.is_disp_ml)
                J_s.yx[k][j][I0.index] -=
                        ml.kappa.x[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][3] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][3]));
            }
          }
          if (I1.apply) {//Perform across I1

            if (!params.is_multilayer) array_ind = I1.index;
            else
              array_ind = (I_tot + 1) * k + I1.index;

            if (k < (K1.index) || params.dimension == TM) {
              E_s.zx[k][j][I1.index] =
                      E_s.zx[k][j][I1.index] +
                      C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
              if (is_cond)
                J_c.zx[k][j][I1.index] -=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
              if (params.is_disp_ml)
                J_s.zx[k][j][I1.index] -=
                        ml.kappa.x[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][6] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][6]));
            }
            if (j < (J1.index)) {
              E_s.yx[k][j][I1.index] =
                      E_s.yx[k][j][I1.index] -
                      C.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
              if (is_cond)
                J_c.yx[k][j][I1.index] +=
                        rho_cond.x[array_ind] * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
              if (params.is_disp_ml)
                J_s.yx[k][j][I1.index] +=
                        ml.kappa.x[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.x[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Isource.real[k - (K0.index)][j - (J0.index)][7] +
                              I * Isource.imag[k - (K0.index)][j - (J0.index)][7]));
            }
          }
        }

      for (k = (K0.index); k <= (K1.index); k++)
        for (i = (I0.index); i <= (I1.index); i++) {
          if (J0.apply) {//Perform across J0
            if (k < (K1.index) || params.dimension == TM) {

              if (!params.is_multilayer) array_ind = J0.index;
              else
                array_ind = (J_tot + 1) * k + J0.index;

              E_s.zy[k][(J0.index)][i] =
                      E_s.zy[k][(J0.index)][i] +
                      C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
              if (is_cond)
                J_c.zy[k][(J0.index)][i] -=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
              if (params.is_disp_ml)
                J_s.zy[k][(J0.index)][i] -=
                        ml.kappa.y[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][2] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][2]));
            }
            if (i < (I1.index)) {
              E_s.xy[k][(J0.index)][i] =
                      E_s.xy[k][(J0.index)][i] -
                      C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
              if (is_cond)
                J_c.xy[k][(J0.index)][i] +=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
              if (params.is_disp_ml)
                J_s.xy[k][(J0.index)][i] +=
                        ml.kappa.y[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][3] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][3]));
            }
          }
          if (J1.apply) {//Perform across J1

            if (!params.is_multilayer) array_ind = J1.index;
            else
              array_ind = (J_tot + 1) * k + J1.index;

            if (k < (K1.index) || params.dimension == TM) {
              E_s.zy[k][(J1.index)][i] =
                      E_s.zy[k][(J1.index)][i] -
                      C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
              if (is_cond)
                J_c.zy[k][(J1.index)][i] +=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
              if (params.is_disp_ml)
                J_s.zy[k][(J1.index)][i] -=
                        ml.kappa.y[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][6] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][6]));
            }
            if (i < (I1.index)) {
              E_s.xy[k][(J1.index)][i] =
                      E_s.xy[k][(J1.index)][i] +
                      C.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
              if (is_cond)
                J_c.xy[k][(J1.index)][i] -=
                        rho_cond.y[array_ind] * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
              if (params.is_disp_ml)
                J_s.xy[k][(J1.index)][i] +=
                        ml.kappa.y[array_ind] * ml.gamma[k] / (2. * params.dt) * C.b.y[array_ind] *
                        real(commonAmplitude * commonPhase *
                             (Jsource.real[k - (K0.index)][i - (I0.index)][7] +
                              I * Jsource.imag[k - (K0.index)][i - (I0.index)][7]));
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
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
              if (is_cond)
                J_c.yz[(K0.index)][j][i] +=
                        rho_cond.z[(K0.index)] * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
              if (params.is_disp_ml)
                J_s.yz[(K0.index)][j][i] -=
                        ml.kappa.z[(K0.index)] * ml.gamma[k] / (2. * params.dt) * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][2]));
            }
            if (i < (I1.index)) {
              E_s.xz[(K0.index)][j][i] =
                      E_s.xz[(K0.index)][j][i] +
                      C.b.z[K0.index] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
              if (is_cond)
                J_c.xz[(K0.index)][j][i] -=
                        rho_cond.z[(K0.index)] * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
              if (params.is_disp_ml)
                J_s.xz[(K0.index)][j][i] +=
                        ml.kappa.z[(K0.index)] * ml.gamma[k] / (2. * params.dt) * C.b.z[K0.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][3]));
            }
          }
          if (K1.apply) {//Perform across K1
            if (j < (J1.index)) {
              E_s.yz[(K1.index)][j][i] =
                      E_s.yz[(K1.index)][j][i] +
                      C.b.z[K1.index] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
              if (is_cond)
                J_c.yz[(K1.index)][j][i] -=
                        rho_cond.z[(K1.index)] * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
              if (params.is_disp_ml)
                J_s.yz[(K1.index)][j][i] +=
                        ml.kappa.z[(K1.index)] * ml.gamma[k] / (2. * params.dt) * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][6] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][6]));
            }
            if (i < (I1.index)) {
              E_s.xz[(K1.index)][j][i] =
                      E_s.xz[(K1.index)][j][i] -
                      C.b.z[K1.index] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
              if (is_cond)
                J_c.xz[(K1.index)][j][i] +=
                        rho_cond.z[(K1.index)] * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
              if (params.is_disp_ml)
                J_s.xz[(K1.index)][j][i] -=
                        ml.kappa.z[(K1.index)] * ml.gamma[k] / (2. * params.dt) * C.b.z[K1.index] *
                        real(commonAmplitude * commonPhase *
                             (Ksource.real[j - (J0.index)][i - (I0.index)][7] +
                              I * Ksource.imag[j - (J0.index)][i - (I0.index)][7]));
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
                                I * Ksource.imag[0][i - (I0.index)][2]) *
                               (-1.0 * I) *
                               exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
                          exp(-1.0 * dcpi *
                              pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
          //Eyz[(int)K0[0]][j][i] = Eyz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + I*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
          if (is_cond)
            J_c.yz[K0.index][j][i] +=
                    rho_cond.z[K0.index] * C.b.z[K0.index] *
                    real((Ksource.real[0][i - (I0.index)][2] +
                          I * Ksource.imag[0][i - (I0.index)][2]) *
                         (-1.0 * I) * exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
                    exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
          //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + I*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
          if (params.is_disp_ml) {
            J_s.yz[K0.index][j][i] -=
                    ml.kappa.z[K0.index] * ml.gamma[K0.index] / (2. * params.dt) *
                    C.b.z[K0.index] *
                    real((Ksource.real[0][i - (I0.index)][2] +
                          I * Ksource.imag[0][i - (I0.index)][2]) *
                         (-1.0 * I) * exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
                    exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
            //Jyz[(int)K0[0]][j][i] -= ml.kappa.z[(int)K0[0]]*ml.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[0][i-((int)I0[0])][2] + I*Ksource.imag[0][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
          }
        }
      } else
        for (j = 0; j < J_tot; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            /*
        if(i==41 & j==41)
        fprintf(stderr,"C.b.z = %.10e, Re(K) = %.10e, Im(K) = %.10e, time_H= %.10e, params.to_l=%.10e, dz/light_v/2=%.10e, hwhm = %.10e, dE=%.10e\n",C.b.z[(int)K0[0]],Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2],Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2],time_H,params.to_l,dz/light_v/2,params.hwhm,C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l + dz/light_v/2.)/(params.hwhm),2)));
      */
            E_s.yz[K0.index][j][i] =
                    E_s.yz[K0.index][j][i] -
                    C.b.z[K0.index] *
                            real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                                  I * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                                 (-1.0 * I) *
                                 exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
                            exp(-1.0 * dcpi *
                                pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
            //Eyz[(int)K0[0]][j][i] = Eyz[(int)K0[0]][j][i] - C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
            if (is_cond)
              J_c.yz[K0.index][j][i] +=
                      rho_cond.z[K0.index] * C.b.z[K0.index] *
                      real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                            I * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                           (-1.0 * I) *
                           exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
                      exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
            //J_c.yz[(int)K0[0]][j][i] += rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
            if (params.is_disp_ml) {
              J_s.yz[K0.index][j][i] -=
                      ml.kappa.z[K0.index] * ml.gamma[K0.index] / (2. * params.dt) *
                      C.b.z[K0.index] *
                      real((Ksource.real[j - (J0.index)][i - (I0.index)][2] +
                            I * Ksource.imag[j - (J0.index)][i - (I0.index)][2]) *
                           (-1.0 * I) *
                           exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
                      exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
              //Jyz[(int)K0[0]][j][i] -= ml.kappa.z[(int)K0[0]]*ml.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][2] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][2])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
            }
          }
      for (j = 0; j < (J_tot + 1); j++)
        for (i = 0; i < I_tot; i++) {
          E_s.xz[K0.index][j][i] =
                  E_s.xz[K0.index][j][i] +
                  C.b.z[K0.index] *
                          real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                                I * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                               (-1.0 * I) *
                               exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2 * dcpi))) *
                          exp(-1.0 * dcpi *
                              pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
          //Exz[(int)K0[0]][j][i] = Exz[(int)K0[0]][j][i] + C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2 ));
          if (is_cond)
            J_c.xz[K0.index][j][i] -=
                    rho_cond.z[K0.index] * C.b.z[K0.index] *
                    real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                          I * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                         (-1.0 * I) * exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2 * dcpi))) *
                    exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
          //J_c.xz[(int)K0[0]][j][i] -= rho_cond.z[(int)K0[0]]*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2 ));
          if (params.is_disp_ml)
            J_s.xz[K0.index][j][i] +=
                    ml.kappa.z[K0.index] * ml.gamma[K0.index] / (2. * params.dt) *
                    C.b.z[K0.index] *
                    real((Ksource.real[j - (J0.index)][i - (I0.index)][3] +
                          I * Ksource.imag[j - (J0.index)][i - (I0.index)][3]) *
                         (-1.0 * I) * exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2 * dcpi))) *
                    exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
          //Jxz[(int)K0[0]][j][i] += ml.kappa.z[(int)K0[0]]*ml.gamma[(int)K0[0]]/(2.*params.dt)*C.b.z[(int)K0[0]]*real((Ksource.real[j-((int)J0[0])][i-((int)I0[0])][3] + I*Ksource.imag[j-((int)J0[0])][i-((int)I0[0])][3])*(-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2 ));
        }
      //fth = real((-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
      H.ft = real((-1.0 * I) * exp(-I * fmod(params.omega_an * (time_H - params.to_l), 2. * dcpi))) *
             exp(-1.0 * dcpi * pow((time_H - params.to_l + dz / light_v / 2.) / (params.hwhm), 2));
      //fth = real((-1.0*I)*exp(-I*fmod(params.omega_an*(time_H - params.to_l),2.*dcpi)))*exp( -1.0*dcpi*pow((time_H - params.to_l)/(params.hwhm),2));
    }
    //fprintf(stderr,"Pos 10:\n");

    //end of source terms
    if (TIME_EXEC) { timer.click(); }

    /********************/
    //begin parallel
#pragma omp parallel default(shared) private(i, j, k, n, k_loc,                                    \
                                             array_ind)//,ca_vec,cb_vec,cc_vec,eh_vec)
    {
      n = omp_get_thread_num();

      if (params.dimension == THREE || params.dimension == TE) {
#ifdef FDFLAG// Use central difference derivatives
//FDTD, Hxz
#pragma omp for
        //Hxz updates
        for (k = 0; k < K_tot; k++)
          for (j = 0; j < J_tot_bound; j++)
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

              if (!materials[k][j][i])
                Hxz[k][j][i] = D.a.z[k_loc] * Hxz[k][j][i] +
                               D.b.z[k_loc] * (Eyx[k + 1][j][i] + Eyz[k + 1][j][i] - Eyx[k][j][i] -
                                             Eyz[k][j][i]);
              else
                Hxz[k][j][i] =
                        Dmaterial.Da.z[materials[k][j][i] - 1] * Hxz[k][j][i] +
                        Dmaterial.Db.z[materials[k][j][i] - 1] *
                                (Eyx[k + 1][j][i] + Eyz[k + 1][j][i] - Eyx[k][j][i] - Eyz[k][j][i]);
            }
//FDTD, Hxz
#else//PSTD, Hxz
#pragma omp for
        //Hxz updates
        for (j = 0; j < J_tot_bound; j++)
          for (i = 0; i < (I_tot + 1); i++) {
            for (k = 0; k < K_tot; k++) {
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

              if (!materials[k][j][i]) {
                ca_vec[n][k] = D.a.z[k_loc];
                cb_vec[n][k] = D.b.z[k_loc];
                //Hxz[k][j][i] = D.a.z[k_loc]*Hxz[k][j][i]+D.b.z[k_loc]*(Eyx[k+1][j][i] + Eyz[k+1][j][i] - Eyx[k][j][i] - Eyz[k][j][i]);
              } else {
                ca_vec[n][k] = Dmaterial.a.z[materials[k][j][i] - 1];
                cb_vec[n][k] = Dmaterial.b.z[materials[k][j][i] - 1];
                //Hxz[k][j][i] = Dmaterial.Da.z[materials[k][j][i]-1]*Hxz[k][j][i]+Dmaterial.Db.z[materials[k][j][i]-1]*(Eyx[k+1][j][i] + Eyz[k+1][j][i] - Eyx[k][j][i] - Eyz[k][j][i]);
              }

              eh_vec[n][k][0] = E_s.yx[k][j][i] + E_s.yz[k][j][i];
              eh_vec[n][k][1] = 0.;
            }
            k = K_tot;
            eh_vec[n][k][0] = E_s.yx[k][j][i] + E_s.yz[k][j][i];
            eh_vec[n][k][1] = 0.;

            first_derivative(eh_vec[n], eh_vec[n], dk_h_z,
                             N_h_z, H_s.xz.plan_f[n], H_s.xz.plan_b[n]);

            for (k = 0; k < K_tot; k++) {
              H_s.xz[k][j][i] = ca_vec[n][k] * H_s.xz[k][j][i] +
                                cb_vec[n][k] *
                                        eh_vec[n][k][0] / ((double) N_h_z);
            }
          }

          //PSTD, Hxz
#endif

#ifdef FDFLAG// Use central difference derivatives
//FDTD, Hxy
#pragma omp for
        //Hxy updates
        for (k = 0; k < K_tot; k++)
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
                Hxy[k][j][i] = D.a.y[array_ind] * Hxy[k][j][i] +
                               D.b.y[array_ind] * (Ezy[k][j][i] + Ezx[k][j][i] - Ezy[k][j + 1][i] -
                                                 Ezx[k][j + 1][i]);
              else
                Hxy[k][j][i] =
                        Dmaterial.Da.y[materials[k][j][i] - 1] * Hxy[k][j][i] +
                        Dmaterial.Db.y[materials[k][j][i] - 1] *
                                (Ezy[k][j][i] + Ezx[k][j][i] - Ezy[k][j + 1][i] - Ezx[k][j + 1][i]);
            }
//FDTD, Hxy
#else//PSTD, Hxy
#pragma omp for
        //Hxy updates
        for (k = 0; k < K_tot; k++)
          for (i = 0; i < (I_tot + 1); i++) {
            for (j = 0; j < J_tot; j++) {
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
              if (!materials[k][j][i]) {
                ca_vec[n][j] = D.a.y[array_ind];
                cb_vec[n][j] = D.b.y[array_ind];
                //		Hxy[k][j][i] = D.a.y[array_ind]*Hxy[k][j][i]+D.b.y[array_ind]*(Ezy[k][j][i] + Ezx[k][j][i] - Ezy[k][j+1][i] - Ezx[k][j+1][i]);
              } else {
                ca_vec[n][j] = Dmaterial.a.y[materials[k][j][i] - 1];
                cb_vec[n][j] = Dmaterial.b.y[materials[k][j][i] - 1];
                //		Hxy[k][j][i] = Dmaterial.Da.y[materials[k][j][i]-1]*Hxy[k][j][i]+Dmaterial.Db.y[materials[k][j][i]-1]*(Ezy[k][j][i] + Ezx[k][j][i] - Ezy[k][j+1][i] - Ezx[k][j+1][i]);
              }

              eh_vec[n][j][0] = E_s.zy[k][j][i] + E_s.zx[k][j][i];
              eh_vec[n][j][1] = 0.;
            }
            j = J_tot;
            eh_vec[n][j][0] = E_s.zy[k][j][i] + E_s.zx[k][j][i];
            eh_vec[n][j][1] = 0.;

            first_derivative(eh_vec[n], eh_vec[n], dk_h_y,
                             N_h_y, H_s.xy.plan_f[n], H_s.xy.plan_b[n]);

            for (j = 0; j < J_tot; j++) {
              H_s.xy[k][j][i] = ca_vec[n][j] * H_s.xy[k][j][i] -
                                cb_vec[n][j] *
                                        eh_vec[n][j][0] / ((double) N_h_y);
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
//PSTD, Hxy
#endif

#ifdef FDFLAG// Use central difference derivatives
//FDTD, Hyx
#pragma omp for
        //Hyx updates
        for (k = 0; k < K_tot; k++)
          for (j = 0; j < J_tot_p1_bound; j++)
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
                Hyx[k][j][i] = D.a.x[array_ind] * Hyx[k][j][i] +
                               D.b.x[array_ind] * (Ezx[k][j][i + 1] + Ezy[k][j][i + 1] -
                                                 Ezx[k][j][i] - Ezy[k][j][i]);
              else {
                Hyx[k][j][i] =
                        Dmaterial.Da.x[materials[k][j][i] - 1] * Hyx[k][j][i] +
                        Dmaterial.Db.x[materials[k][j][i] - 1] *
                                (Ezx[k][j][i + 1] + Ezy[k][j][i + 1] - Ezx[k][j][i] - Ezy[k][j][i]);
              }
            }
//FDTD, Hyx
#else//PSTD, Hyx
#pragma omp for
        //Hyx updates
        for (k = 0; k < K_tot; k++)
          for (j = 0; j < J_tot_p1_bound; j++) {
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
                ca_vec[n][i] = D.a.x[array_ind];
                cb_vec[n][i] = D.b.x[array_ind];
                //		Hyx[k][j][i] = D.a.x[array_ind]*Hyx[k][j][i]+D.b.x[array_ind]*(Ezx[k][j][i+1] + Ezy[k][j][i+1] - Ezx[k][j][i] - Ezy[k][j][i]);
              } else {
                ca_vec[n][i] = Dmaterial.a.x[materials[k][j][i] - 1];
                cb_vec[n][i] = Dmaterial.b.x[materials[k][j][i] - 1];
                //	Hyx[k][j][i] = Dmaterial.Da.x[materials[k][j][i]-1]*Hyx[k][j][i]+Dmaterial.Db.x[materials[k][j][i]-1]*(Ezx[k][j][i+1] + Ezy[k][j][i+1] - Ezx[k][j][i] - Ezy[k][j][i]);
              }

              eh_vec[n][i][0] = E_s.zx[k][j][i] + E_s.zy[k][j][i];
              eh_vec[n][i][1] = 0.;
            }
            i = I_tot;
            eh_vec[n][i][0] = E_s.zx[k][j][i] + E_s.zy[k][j][i];
            eh_vec[n][i][1] = 0.;

            first_derivative(eh_vec[n], eh_vec[n], dk_h_x,
                             N_h_x, H_s.yx.plan_f[n], H_s.yx.plan_b[n]);

            for (i = 0; i < I_tot; i++) {
              H_s.yx[k][j][i] = ca_vec[n][i] * H_s.yx[k][j][i] +
                                cb_vec[n][i] *
                                        eh_vec[n][i][0] / ((double) N_h_x);
            }
          }
//PSTD, Hyx
#endif

#ifdef FDFLAG// Use central difference derivatives
//FDTD, Hyz
#pragma omp for
        //Hyz updates
        for (k = 0; k < K_tot; k++) {
          for (j = 0; j < J_tot_p1_bound; j++)
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
              if (!materials[k][j][i]) {
                /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                Hyz[k][j][i] = D.a.z[k_loc] * Hyz[k][j][i] +
                               D.b.z[k_loc] * (Exy[k][j][i] + Exz[k][j][i] - Exy[k + 1][j][i] -
                                             Exz[k + 1][j][i]);
              } else {
                /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                Hyz[k][j][i] =
                        Dmaterial.Da.z[materials[k][j][i] - 1] * Hyz[k][j][i] +
                        Dmaterial.Db.z[materials[k][j][i] - 1] *
                                (Exy[k][j][i] + Exz[k][j][i] - Exy[k + 1][j][i] - Exz[k + 1][j][i]);
              }
            }
        }
//FDTD, Hyz
#else//PSTD, Hyz
        //#pragma omp for
        //Hyz updates
        for (j = 0; j < J_tot_p1_bound; j++)
#pragma omp for
          for (i = 0; i < I_tot; i++) {
            for (k = 0; k < K_tot; k++) {
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
              if (!materials[k][j][i]) {
                ca_vec[n][k] = D.a.z[k_loc];
                cb_vec[n][k] = D.b.z[k_loc];
                /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,D.a.z[k_loc], D.b.z[k_loc]);*/
                //Hyz[k][j][i] = D.a.z[k_loc]*Hyz[k][j][i]+D.b.z[k_loc]*(Exy[k][j][i] + Exz[k][j][i] - Exy[k+1][j][i] - Exz[k+1][j][i]);
              } else {
                ca_vec[n][k] = Dmaterial.a.z[materials[k][j][i] - 1];
                cb_vec[n][k] = Dmaterial.b.z[materials[k][j][i] - 1];
                /*if(tind==0)
        fprintf(stdout,"%d %d %e %e\n",i,k,Dmaterial.Da.z[materials[k][j][i]-1],Dmaterial.Db.z[materials[k][j][i]-1]);*/
                //Hyz[k][j][i] = Dmaterial.Da.z[materials[k][j][i]-1]*Hyz[k][j][i]+Dmaterial.Db.z[materials[k][j][i]-1]*(Exy[k][j][i] + Exz[k][j][i] - Exy[k+1][j][i] - Exz[k+1][j][i]);
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

            first_derivative(eh_vec[n], eh_vec[n], dk_h_z,
                             N_h_z, H_s.yz.plan_f[n], H_s.yz.plan_b[n]);

            for (k = 0; k < K_tot; k++) {
              H_s.yz[k][j][i] = ca_vec[n][k] * H_s.yz[k][j][i] -
                                cb_vec[n][k] *
                                        eh_vec[n][k][0] / ((double) N_h_z);
            }
          }
//PSTD, Hyz
#endif
      }//(params.dimension==THREE || params.dimension==TE)
      else {

#pragma omp for
        for (k = 0; k <= K_tot; k++)
          for (j = 0; j < J_tot; j++)
            for (i = 0; i < (I_tot + 1); i++)
              if (!materials[k][j][i]) H_s.xz[k][j][i] = 0.;
              else
                H_s.xz[k][j][i] = 0.;

#pragma omp for
        //Hxy update
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
        //Hyx update
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

      if (params.dimension == THREE || params.dimension == TE) {
#ifdef FDFLAG// Use central difference derivatives
//FDTD, Hzy
#pragma omp for
        //Hzy update
        for (k = 0; k < (K_tot + 1); k++)
          for (j = 0; j < J_tot; j++)
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
              if (!params.is_multilayer) array_ind = j;
              else
                array_ind = (J_tot + 1) * k_loc + j;
              if (!materials[k][j][i])
                Hzy[k][j][i] = D.a.y[array_ind] * Hzy[k][j][i] +
                               D.b.y[array_ind] * (Exy[k][j + 1][i] + Exz[k][j + 1][i] -
                                                 Exy[k][j][i] - Exz[k][j][i]);
              else
                Hzy[k][j][i] =
                        Dmaterial.Da.y[materials[k][j][i] - 1] * Hzy[k][j][i] +
                        Dmaterial.Db.y[materials[k][j][i] - 1] *
                                (Exy[k][j + 1][i] + Exz[k][j + 1][i] - Exy[k][j][i] - Exz[k][j][i]);
            }
//FDTD, Hzy
#else//PSTD, Hzy
#pragma omp for
        //Hzy update
        for (k = 0; k < (K_tot + 1); k++)
          for (i = 0; i < I_tot; i++) {
            for (j = 0; j < J_tot; j++) {
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
              if (!materials[k][j][i]) {
                ca_vec[n][j] = D.a.y[array_ind];
                cb_vec[n][j] = D.b.y[array_ind];
                //	      Hzy[k][j][i] = D.a.y[array_ind]*Hzy[k][j][i]+D.b.y[array_ind]*(Exy[k][j+1][i] + Exz[k][j+1][i] - Exy[k][j][i] - Exz[k][j][i]);
              } else {
                ca_vec[n][j] = Dmaterial.a.y[materials[k][j][i] - 1];
                cb_vec[n][j] = Dmaterial.b.y[materials[k][j][i] - 1];
                //	      Hzy[k][j][i] = Dmaterial.Da.y[materials[k][j][i]-1]*Hzy[k][j][i]+Dmaterial.Db.y[materials[k][j][i]-1]*(Exy[k][j+1][i] + Exz[k][j+1][i] - Exy[k][j][i] - Exz[k][j][i]);
              }

              eh_vec[n][j][0] = E_s.xy[k][j][i] + E_s.xz[k][j][i];
              eh_vec[n][j][1] = 0.;
            }
            j = J_tot;
            eh_vec[n][j][0] = E_s.xy[k][j][i] + E_s.xz[k][j][i];
            eh_vec[n][j][1] = 0.;

            first_derivative(eh_vec[n], eh_vec[n], dk_h_y, N_h_y, H_s.zy.plan_f[n], H_s.zy.plan_b[n]);

            for (j = 0; j < J_tot; j++) {
              H_s.zy[k][j][i] = ca_vec[n][j] * H_s.zy[k][j][i] +
                                cb_vec[n][j] *
                                        eh_vec[n][j][0] / ((double) N_h_y);
            }
          }
          //PSTD, Hzy
#endif


#ifdef FDFLAG// Use central difference derivatives
//FDTD, Hzx
#pragma omp for
        //Hzx update
        for (k = 0; k < (K_tot + 1); k++)
          for (j = 0; j < J_tot_bound; j++)
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
                Hzx[k][j][i] = D.a.x[array_ind] * Hzx[k][j][i] +
                               D.b.x[array_ind] * (Eyx[k][j][i] + Eyz[k][j][i] - Eyx[k][j][i + 1] -
                                                 Eyz[k][j][i + 1]);
              else
                Hzx[k][j][i] =
                        Dmaterial.Da.x[materials[k][j][i] - 1] * Hzx[k][j][i] +
                        Dmaterial.Db.x[materials[k][j][i] - 1] *
                                (Eyx[k][j][i] + Eyz[k][j][i] - Eyx[k][j][i + 1] - Eyz[k][j][i + 1]);
            }
//FDTD, Hzx
#else//PSTD, Hzx
#pragma omp for
        //Hzx update
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
                //		Hzx[k][j][i] = D.a.x[array_ind]*Hzx[k][j][i]+D.b.x[array_ind]*(Eyx[k][j][i] + Eyz[k][j][i] - Eyx[k][j][i+1] - Eyz[k][j][i+1]);
                ca_vec[n][i] = D.a.x[array_ind];
                cb_vec[n][i] = D.b.x[array_ind];
              } else {
                //		Hzx[k][j][i] = Dmaterial.Da.x[materials[k][j][i]-1]*Hzx[k][j][i]+Dmaterial.Db.x[materials[k][j][i]-1]*(Eyx[k][j][i] + Eyz[k][j][i] - Eyx[k][j][i+1] - Eyz[k][j][i+1]);
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
//PSTD, Hzx
#endif
      }//(params.dimension==THREE || params.dimension==TE)
    }  //end parallel
    if (TIME_EXEC) { timer.click(); }

    //fprintf(stderr,"Pos 11b:\n");
    //update terms for self consistency across scattered/total interface - E updates
    if (params.source_mode == SourceMode::steadystate) {//steadystate
      complex<double> commonPhase = exp(-I * fmod(params.omega_an * time_E, 2. * dcpi));
      double commonAmplitude = linearRamp(time_E, 1. / (params.omega_an / (2 * dcpi)), ramp_width);
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
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][0]));
            if (k < (K1.index) || params.dimension == TM)
              H_s.yx[k][j][(I0.index) - 1] =
                      H_s.yx[k][j][(I0.index) - 1] -
                      D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][1] +
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][1]));
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
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][4]));
            if (k < (K1.index) || params.dimension == TM)
              H_s.yx[k][j][(I1.index)] =
                      H_s.yx[k][j][(I1.index)] +
                      D.b.x[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Isource.real[k - (K0.index)][j - (J0.index)][5] +
                                    I * Isource.imag[k - (K0.index)][j - (J0.index)][5]));
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
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][0]));

            if (k < (K1.index) || params.dimension == TM)
              H_s.xy[k][(J0.index) - 1][i] =
                      H_s.xy[k][(J0.index) - 1][i] +
                      D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][1] +
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][1]));
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
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][4]));
            if (k < (K1.index) || params.dimension == TM)
              H_s.xy[k][(J1.index)][i] =
                      H_s.xy[k][(J1.index)][i] -
                      D.b.y[array_ind] *
                              real(commonAmplitude * commonPhase *
                                   (Jsource.real[k - (K0.index)][i - (I0.index)][5] +
                                    I * Jsource.imag[k - (K0.index)][i - (I0.index)][5]));
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
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][0]));
            if (j < (J1.index))
              H_s.xz[(K0.index) - 1][j][i] =
                      H_s.xz[(K0.index) - 1][j][i] -
                      D.b.z[(K0.index) - 1] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][1] +
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][1]));
          }
          if (K1.apply) {//Perform across K1
            if (i < (I1.index))
              H_s.yz[(K1.index)][j][i] =
                      H_s.yz[(K1.index)][j][i] -
                      D.b.z[(K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][4] +
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][4]));
            if (j < (J1.index))
              H_s.xz[(K1.index)][j][i] =
                      H_s.xz[(K1.index)][j][i] +
                      D.b.z[(K1.index)] *
                              real(commonAmplitude * commonPhase *
                                   (Ksource.real[j - (J0.index)][i - (I0.index)][5] +
                                    I * Ksource.imag[j - (J0.index)][i - (I0.index)][5]));
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
                                I * Ksource.imag[0][i - (I0.index)][1]) *
                               (-1. * I) *
                               exp(-I * fmod(params.omega_an * (time_E - params.to_l), 2 * dcpi))) *
                          exp(-1. * dcpi * pow((time_E - params.to_l) / (params.hwhm), 2.));
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
                                I * Ksource.imag[0][i - (I0.index)][0]) *
                               (-1. * I) *
                               exp(-I * fmod(params.omega_an * (time_E - params.to_l), 2 * dcpi))) *
                          exp(-1. * dcpi * pow((time_E - params.to_l) / (params.hwhm), 2.));
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
                                  I * Ksource.imag[j - (J0.index)][i - (I0.index)][1]) *
                                 (-1. * I) *
                                 exp(-I * fmod(params.omega_an * (time_E - params.to_l), 2 * dcpi))) *
                            exp(-1. * dcpi * pow((time_E - params.to_l) / (params.hwhm), 2.));
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
                                  I * Ksource.imag[j - (J0.index)][i - (I0.index)][0]) *
                                 (-1. * I) *
                                 exp(-I * fmod(params.omega_an * (time_E - params.to_l), 2 * dcpi))) *
                            exp(-1. * dcpi * pow((time_E - params.to_l) / (params.hwhm), 2.));
            //broadband source term
            if (params.exi_present)
              H_s.yz[(K0.index) - 1][j][i] =
                      H_s.yz[(K0.index) - 1][j][i] + D.b.z[(K0.index) - 1] * Ei.x[tind][j][i];
          }
        //fprintf(stderr,"Pos 11i\n");
      }
      E.ft = real((-1. * I) * exp(-I * fmod(params.omega_an * (time_E - params.to_l), 2 * dcpi))) *
             exp(-1. * dcpi * pow((time_E - params.to_l) / (params.hwhm), 2.));
      //fprintf(stderr,"Pos 11j\n");
    }
    if (TIME_EXEC) { timer.click(); }

    if (params.exphasorssurface || params.exphasorsvolume || params.exdetintegral || (campssample.n_vertices() > 0)) {
      if (params.source_mode == SourceMode::steadystate) {
        E.add_to_angular_norm(tind, Nsteps, params);
        H.add_to_angular_norm(tind, Nsteps, params);

        for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
          extractPhasorENorm(&E_norm[ifx], E.ft, tind, f_ex_vec[ifx] * 2 * dcpi, params.dt, Nsteps);
          extractPhasorHNorm(&H_norm[ifx], H.ft, tind, f_ex_vec[ifx] * 2 * dcpi, params.dt, Nsteps);
        }
      } else {
        if ((tind - params.start_tind) % params.Np == 0) {

          E.add_to_angular_norm(tind, params.Npe, params);
          H.add_to_angular_norm(tind, params.Npe, params);

          for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
            extractPhasorENorm(&E_norm[ifx], E.ft, tind, f_ex_vec[ifx] * 2 * dcpi, params.dt, params.Npe);
            extractPhasorHNorm(&H_norm[ifx], H.ft, tind, f_ex_vec[ifx] * 2 * dcpi, params.dt, params.Npe);
          }
        }
      }
    }
    if (TIME_EXEC) { timer.click(); }

    if ((((double) time(NULL)) - t0) > 1) {

      maxfield = 0.;
      for (k = 0; k < (K_tot + 1); k++) {
        for (j = 0; j < (J_tot + 1); j++) {
          for (i = 0; i < (I_tot + 1); i++) {
            tempfield = fabs(E_s.xy[k][j][i] + E_s.xz[k][j][i]);
            if (maxfield < tempfield) { maxfield = tempfield; }
            tempfield = fabs(E_s.yx[k][j][i] + E_s.yz[k][j][i]);
            if (maxfield < tempfield) { maxfield = tempfield; }
            tempfield = fabs(E_s.zx[k][j][i] + E_s.zy[k][j][i]);
            if (maxfield < tempfield) { maxfield = tempfield; }
            tempfield = fabs(H_s.xy[k][j][i] + H_s.xz[k][j][i]);
            if (maxfield < tempfield) { maxfield = tempfield; }
            tempfield = fabs(H_s.yx[k][j][i] + H_s.yz[k][j][i]);
            if (maxfield < tempfield) { maxfield = tempfield; }
            tempfield = fabs(H_s.zx[k][j][i] + H_s.zy[k][j][i]);
            if (maxfield < tempfield) { maxfield = tempfield; }
          }
        }
      }

      fprintf(stdout, "Iterating: %d %e\n", tind, maxfield);
      //fprintf(stderr,"Post-iter 1\n");
      //     fprintf(stdout,"Iterating: %d\n",tind);
      t0 = double(time(NULL));
      //fprintf(stderr,"Post-iter 2\n");
    }
    //fprintf(stderr,"Post-iter 3\n");
    if ((params.source_mode == SourceMode::steadystate) && (tind == (params.Nt - 1)) && (params.run_mode == RunMode::complete) &&
        params.exphasorsvolume) {
      fprintf(stdout, "Iteration limit reached, setting output fields to last complete DFT\n");
      copyPhasors(E_copy, E, (int) mxGetNumberOfElements((mxArray *) plhs[0]));
    }
    //fprintf(stderr,"Post-iter 4\n");
    fflush(stdout);
    //fprintf(stderr,"Post-iter 5\n");
    //fprintf(stderr,"%s %d %d\n",tdfdirstr, strcmp(tdfdirstr,""),are_equal(tdfdirstr,""));
    if (params.has_tdfdir && (tind % params.Np) == 0) {
      fprintf(stderr,"Saving field\n");
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
       matPutVariable(toutfile, "fdtdgrid", (mxArray *)prhs[0]);
       matClose(toutfile);
       }
    */
    /*write out fdtdgrid to a file*/

  }//end of main iteration loop
  if (TIME_MAIN_LOOP) {
    //fprintf(stderr,"Post-iter 7\n");
    main_loop_timer.end();
    //fprintf(stderr,"Post-iter 8\n");
    fprintf(stdout, "# Time elasped in main loop: %e\n", main_loop_timer.delta_seconds());
    //fprintf(stderr,"Post-iter 9\n");
  }
  //save state of fdtdgrid

  //fprintf(stderr,"Pos 12\n");
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    E.normalise_volume();
    H.normalise_volume();
  }

  //fprintf(stderr,"Pos 13\n");
  if (params.run_mode == RunMode::complete && params.exphasorssurface)
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      normaliseSurface(surface_EHr[ifx], surface_EHi[ifx], surface_vertices, n_surface_vertices,
                       E_norm[ifx], H_norm[ifx]);
      //fprintf(stderr,"E_norm[%d]: %e %e\n",ifx,real(E_norm[ifx]),imag(E_norm[ifx]));
    }

  if (params.run_mode == RunMode::complete && (campssample.n_vertices() > 0))
    for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
      normaliseVertices(camplitudesR[ifx], camplitudesI[ifx], campssample, E_norm[ifx], H_norm[ifx]);
      fprintf(stderr, "E_norm[%d]: %e %e\n", ifx, real(E_norm[ifx]), imag(E_norm[ifx]));
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
  maxfield = 0.0;
  for (k = 0; k < (K_tot + 1); k++) {
    for (j = 0; j < (J_tot + 1); j++) {
      for (i = 0; i < (I_tot + 1); i++) {
        tempfield = fabs(E_s.xy[k][j][i] + E_s.xz[k][j][i]);
        if (maxfield < tempfield) { maxfield = tempfield; }
        tempfield = fabs(E_s.yx[k][j][i] + E_s.yz[k][j][i]);
        if (maxfield < tempfield) { maxfield = tempfield; }
        tempfield = fabs(E_s.zx[k][j][i] + E_s.zy[k][j][i]);
        if (maxfield < tempfield) { maxfield = tempfield; }
        tempfield = fabs(H_s.xy[k][j][i] + H_s.xz[k][j][i]);
        if (maxfield < tempfield) { maxfield = tempfield; }
        tempfield = fabs(H_s.yx[k][j][i] + H_s.yz[k][j][i]);
        if (maxfield < tempfield) { maxfield = tempfield; }
        tempfield = fabs(H_s.zx[k][j][i] + H_s.zy[k][j][i]);
        if (maxfield < tempfield) { maxfield = tempfield; }
      }
    }
  }

  //fprintf(stderr,"Pos 15\n");
  //noe set the output
  ndims = 2;
  dims[0] = 1;
  dims[1] = 1;
  plhs[25] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  *mxGetPr((mxArray *) plhs[25]) = maxfield;

  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    setGridLabels(input_grid_labels, output_grid_labels, E.il, E.iu, E.jl, E.ju, E.kl, E.ku);
  }
  
  auto interp_output_grid_labels = GridLabels();

  //fprintf(stderr,"Pos 15_m1\n");
  if (params.run_mode == RunMode::complete && params.exphasorsvolume) {
    //now interpolate over the extracted phasors
    if (params.dimension == THREE) {
      fprintf(stderr, "mxInterpolateFieldCentralE: %d %d %d \n", E.I_tot - 2,
              E.J_tot - 2, E.K_tot - 2);
      //fprintf(stderr,"Pos 15_m1a\n");
      mxInterpolateFieldCentralE(plhs[0], plhs[1], plhs[2], &plhs[13], &plhs[14], &plhs[15], 2,
                                 E.I_tot - 2, 2, E.J_tot - 2, 2,
                                 E.K_tot - 2);
      //fprintf(stderr,"Pos 15_m1b\n");

    } else if (params.dimension == TE)
      mxInterpolateFieldCentralE_TE(plhs[0], plhs[1], plhs[2], &plhs[13], &plhs[14], &plhs[15], 2,
                                    E.I_tot - 2, 2, E.J_tot - 2, 0, 0);
    else
      mxInterpolateFieldCentralE_TM(plhs[0], plhs[1], plhs[2], &plhs[13], &plhs[14], &plhs[15], 2,
                                    E.I_tot - 2, 2, E.J_tot - 2, 0, 0);
    if (params.dimension == THREE)
      mxInterpolateFieldCentralH(plhs[3], plhs[4], plhs[5], &plhs[16], &plhs[17], &plhs[18], 2,
                                 E.I_tot - 2, 2, E.J_tot - 2, 2,
                                 E.K_tot - 2);
    else if (params.dimension == TE)
      mxInterpolateFieldCentralH_TE(plhs[3], plhs[4], plhs[5], &plhs[16], &plhs[17], &plhs[18], 2,
                                    E.I_tot - 2, 2, E.J_tot - 2, 0, 0);
    else
      mxInterpolateFieldCentralH_TM(plhs[3], plhs[4], plhs[5], &plhs[16], &plhs[17], &plhs[18], 2,
                                    E.I_tot - 2, 2, E.J_tot - 2, 0, 0);

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
      //fprintf(stderr,"Pos 15a-1\n");
      setGridLabels(output_grid_labels, interp_output_grid_labels, 2, E.I_tot - 2, 2,
                    E.J_tot - 2, 2, E.K_tot - 2);
    } else
      setGridLabels(output_grid_labels, interp_output_grid_labels, 2, E.I_tot - 2, 2,
                    E.J_tot - 2, 0, 0);
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
                                   params.phasorinc, &dummy_vertex_list,
                                   &mx_surface_facets);
    mxDestroyArray(dummy_vertex_list);
    mxArray *vertex_list;
    double **vertex_list_ptr;
    ndims = 2;
    dims[0] = n_surface_vertices;
    dims[1] = 3;
    vertex_list = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
    vertex_list_ptr = cast_matlab_2D_array(mxGetPr((mxArray *) vertex_list), dims[0], dims[1]);

    //now populate the vertex list
    for (i = 0; i < n_surface_vertices; i++) {

      vertex_list_ptr[0][i] = input_grid_labels.x[surface_vertices[0][i]];
      vertex_list_ptr[1][i] = input_grid_labels.y[surface_vertices[1][i]];
      vertex_list_ptr[2][i] = input_grid_labels.z[surface_vertices[2][i]];
    }
    //assign outputs
    plhs[22] = vertex_list;
    plhs[23] = mx_surface_amplitudes;
    plhs[24] = mx_surface_facets;

    free_cast_matlab_2D_array(vertex_list_ptr);
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
  /*Free the additional data structures used to cast the matlab arrays*/
  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    free_cast_matlab_2D_array(surface_vertices);
    free_cast_matlab_3D_array(surface_EHr, f_ex_vec.size());
    free_cast_matlab_3D_array(surface_EHi, f_ex_vec.size());

    mxDestroyArray(mx_surface_vertices);
  }

  if (campssample.n_vertices() > 0) {
    free_cast_matlab_3D_array(camplitudesR, f_ex_vec.size());
    free_cast_matlab_3D_array(camplitudesI, f_ex_vec.size());
  }

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
  if (I0.apply || I1.apply) {
    free_cast_matlab_3D_array(Isource.imag, (K1.index - K0.index + 1));
    free_cast_matlab_3D_array(Isource.real, (K1.index - K0.index + 1));
  }
  if (J0.apply || J1.apply) {
    free_cast_matlab_3D_array(Jsource.imag, (K1.index - K0.index + 1));
    free_cast_matlab_3D_array(Jsource.real, (K1.index - K0.index + 1));
  }
  if (K0.apply || K1.apply) {
    free_cast_matlab_3D_array(Ksource.imag, (J1.index - J0.index + 1));
    free_cast_matlab_3D_array(Ksource.real, (J1.index - J0.index + 1));
  }

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

#ifndef FDFLAG// Using PS

  fftw_free(dk_e_x);
  fftw_free(dk_e_y);
  fftw_free(dk_e_z);

  fftw_free(dk_h_x);
  fftw_free(dk_h_y);
  fftw_free(dk_h_z);
#endif

  free(E_norm);
  free(H_norm);
  free(dims);
  free(label_dims);

  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete) {
    mxDestroyArray(dummy_array[0]);
    mxDestroyArray(dummy_array[1]);
    mxDestroyArray(dummy_array[2]);
  }

  //must destroy mx_surface_amplitudes
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

void normaliseSurface(double **surface_EHr, double **surface_EHi, int **surface_vertices,
                      int n_surface_vertices, complex<double> Enorm, complex<double> Hnorm) {
  double norm_r, norm_i, denom, temp_r, temp_i;

  norm_r = real(Enorm);
  norm_i = imag(Enorm);
  denom = norm_r * norm_r + norm_i * norm_i;

  for (int vindex = 0; vindex < n_surface_vertices; vindex++)
    for (int i = 0; i < 3; i++) {
      temp_r = surface_EHr[i][vindex];
      temp_i = surface_EHi[i][vindex];

      surface_EHr[i][vindex] = (norm_r * temp_r + norm_i * temp_i) / denom;
      surface_EHi[i][vindex] = (norm_r * temp_i - norm_i * temp_r) / denom;
    }

  norm_r = real(Hnorm);
  norm_i = imag(Hnorm);
  denom = norm_r * norm_r + norm_i * norm_i;

  for (int vindex = 0; vindex < n_surface_vertices; vindex++)
    for (int i = 3; i < 6; i++) {
      temp_r = surface_EHr[i][vindex];
      temp_i = surface_EHi[i][vindex];

      surface_EHr[i][vindex] = (norm_r * temp_r + norm_i * temp_i) / denom;
      surface_EHi[i][vindex] = (norm_r * temp_i - norm_i * temp_r) / denom;
    }
}

void normaliseVertices(double **EHr, double **EHi, ComplexAmplitudeSample &campssample, complex<double> Enorm, complex<double> Hnorm) {

  for (int i = 0; i < 6; i++) {
    
    auto norm = i < 3 ? Enorm : Hnorm;
    double norm_r = real(norm);
    double norm_i = imag(norm);
    double denom = norm_r * norm_r + norm_i * norm_i;

    auto ii = campssample.components.index(i + 1);
    if (ii >= 0) {
      for (int vindex = 0; vindex < campssample.n_vertices(); vindex++){

        double temp_r = EHr[ii][vindex];
        double temp_i = EHi[ii][vindex];

        EHr[ii][vindex] = (norm_r * temp_r + norm_i * temp_i) / denom;
        EHi[ii][vindex] = (norm_r * temp_i - norm_i * temp_r) / denom;
      }
    }
  }
}

void extractPhasorENorm(complex<double> *Enorm, double ft, int n, double omega, double dt, int Nt) {
  *Enorm += ft * exp(fmod(omega * ((double) (n + 1)) * dt, 2 * dcpi) * I) * 1. / ((double) Nt);
}

void extractPhasorHNorm(complex<double> *Hnorm, double ft, int n, double omega, double dt, int Nt) {
  *Hnorm += ft * exp(fmod(omega * ((double) n + 0.5) * dt, 2 * dcpi) * I) * 1. / ((double) Nt);
}


void extractPhasorsSurface(double **surface_EHr, double **surface_EHi, ElectricSplitField &E,
                           MagneticSplitField &H, int **surface_vertices, int n_surface_vertices,
                           int n, double omega, int Nt, int J_tot, SimulationParameters &params) {
  int vindex;
  double Ex, Ey, Ez, Hx, Hy, Hz;
  complex<double> phaseTermE, phaseTermH, subResultE, subResultH, cphaseTermE, cphaseTermH;

  phaseTermE = fmod(omega * ((double) n) * params.dt, 2 * dcpi);
  phaseTermH = fmod(omega * ((double) n + 0.5) * params.dt, 2 * dcpi);

  cphaseTermH = exp(phaseTermH * I) * 1. / ((double) Nt);
  cphaseTermE = exp(phaseTermE * I) * 1. / ((double) Nt);

  //loop over every vertex in the list
#pragma omp parallel default(shared) private(Ex, Ey, Ez, Hx, Hy, Hz, phaseTermE, phaseTermH,       \
                                             subResultE, subResultH, vindex)
  {
#pragma omp for
    for (vindex = 0; vindex < n_surface_vertices; vindex++) {
      //    fprintf(stderr,"vindex: %d: (%d %d %d)\n",vindex,surface_vertices[0][vindex],surface_vertices[1][vindex],surface_vertices[2][vindex]);
      if (params.dimension == THREE)
        if (J_tot == 0) {
          interpolateTimeDomainFieldCentralE_2Dy(
                  E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, surface_vertices[0][vindex],
                  surface_vertices[1][vindex], surface_vertices[2][vindex], &Ex, &Ey, &Ez);
        } else if (params.interp_method == cubic)
          interpolateTimeDomainFieldCentralE(
                  E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, surface_vertices[0][vindex],
                  surface_vertices[1][vindex], surface_vertices[2][vindex], &Ex, &Ey, &Ez);
        else
          interpolateTimeDomainFieldCentralEBandLimited(
                  E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, surface_vertices[0][vindex],
                  surface_vertices[1][vindex], surface_vertices[2][vindex], &Ex, &Ey, &Ez);
      else if (params.dimension == TE)
        interpolateTimeDomainFieldCentralE_TE(
                E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, surface_vertices[0][vindex],
                surface_vertices[1][vindex], surface_vertices[2][vindex], &Ex, &Ey, &Ez);
      else
        interpolateTimeDomainFieldCentralE_TM(
                E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, surface_vertices[0][vindex],
                surface_vertices[1][vindex], surface_vertices[2][vindex], &Ex, &Ey, &Ez);
      //    fprintf(stderr,"1st interp donezn");
      if (params.dimension == THREE)
        if (J_tot == 0) {
          interpolateTimeDomainFieldCentralH_2Dy(
                  H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, surface_vertices[0][vindex],
                  surface_vertices[1][vindex], surface_vertices[2][vindex], &Hx, &Hy, &Hz);
        } else if (params.interp_method == cubic)
          interpolateTimeDomainFieldCentralH(
                  H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, surface_vertices[0][vindex],
                  surface_vertices[1][vindex], surface_vertices[2][vindex], &Hx, &Hy, &Hz);
        else
          interpolateTimeDomainFieldCentralHBandLimited(
                  H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, surface_vertices[0][vindex],
                  surface_vertices[1][vindex], surface_vertices[2][vindex], &Hx, &Hy, &Hz);
      else if (params.dimension == TE)
        interpolateTimeDomainFieldCentralH_TE(
                H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, surface_vertices[0][vindex],
                surface_vertices[1][vindex], surface_vertices[2][vindex], &Hx, &Hy, &Hz);
      else
        interpolateTimeDomainFieldCentralH_TM(
                H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, surface_vertices[0][vindex],
                surface_vertices[1][vindex], surface_vertices[2][vindex], &Hx, &Hy, &Hz);
      //    fprintf(stderr,"2nd interp donezn");

      /*Ex and Hx*/
      subResultH = Hx * cphaseTermH;//exp(phaseTermH * I) * 1./((double) Nt);
      subResultE = Ex * cphaseTermE;//exp(phaseTermE * I) * 1./((double) Nt);

      //now update the master array
      surface_EHr[0][vindex] = surface_EHr[0][vindex] + real(subResultE);
      surface_EHi[0][vindex] = surface_EHi[0][vindex] + imag(subResultE);

      surface_EHr[3][vindex] = surface_EHr[3][vindex] + real(subResultH);
      surface_EHi[3][vindex] = surface_EHi[3][vindex] + imag(subResultH);

      /*Ey and Hy*/
      subResultH = Hy * cphaseTermH;//exp(phaseTermH * I) * 1./((double) Nt);
      subResultE = Ey * cphaseTermE;//exp(phaseTermE * I) * 1./((double) Nt);

      //now update the master array
      surface_EHr[1][vindex] = surface_EHr[1][vindex] + real(subResultE);
      surface_EHi[1][vindex] = surface_EHi[1][vindex] + imag(subResultE);

      surface_EHr[4][vindex] = surface_EHr[4][vindex] + real(subResultH);
      surface_EHi[4][vindex] = surface_EHi[4][vindex] + imag(subResultH);


      /*Ez and Hz*/
      subResultH = Hz * cphaseTermH;//exp(phaseTermH * I) * 1./((double) Nt);
      subResultE = Ez * cphaseTermE;//exp(phaseTermE * I) * 1./((double) Nt);

      //now update the master array
      surface_EHr[2][vindex] = surface_EHr[2][vindex] + real(subResultE);
      surface_EHi[2][vindex] = surface_EHi[2][vindex] + imag(subResultE);

      surface_EHr[5][vindex] = surface_EHr[5][vindex] + real(subResultH);
      surface_EHi[5][vindex] = surface_EHi[5][vindex] + imag(subResultH);
    }
  }//end parallel region
}

void extractPhasorsVertices(double **EHr, double **EHi, ElectricSplitField &E, MagneticSplitField &H,
                            ComplexAmplitudeSample &campssample, int n, double omega,
                            double dt, int Nt, int dimension, int J_tot, int intmethod) {

  int vindex, i, j, k;
  double Ex, Ey, Ez, Hx, Hy, Hz;
  complex<double> cphaseTermE, cphaseTermH;

  auto phaseTermE = fmod(omega * ((double) n) * dt, 2 * dcpi);
  auto phaseTermH = fmod(omega * ((double) n + 0.5) * dt, 2 * dcpi);

  cphaseTermH = exp(phaseTermH * I) * 1. / ((double) Nt);
  cphaseTermE = exp(phaseTermE * I) * 1. / ((double) Nt);

#pragma omp parallel default(none) \
        shared(E, H, EHr, EHi, campssample) \
        private(Ex, Ey, Ez, Hx, Hy, Hz, vindex, i, j, k) \
        firstprivate(cphaseTermH, cphaseTermE, dimension, J_tot, intmethod)
  {
#pragma omp for
    for (vindex = 0; vindex < campssample.n_vertices(); vindex++) {   // loop over every vertex

      i = campssample.vertices[0][vindex];
      j = campssample.vertices[1][vindex];
      k = campssample.vertices[2][vindex];

      if (dimension == THREE)
        if (J_tot == 0) {
          interpolateTimeDomainFieldCentralE_2Dy(E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, i, j, k, &Ex,
                                                 &Ey, &Ez);
        } else if (intmethod == 1)
          interpolateTimeDomainFieldCentralE(E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, i, j, k, &Ex, &Ey,
                                             &Ez);
        else
          interpolateTimeDomainFieldCentralEBandLimited(E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, i, j, k,
                                                        &Ex, &Ey, &Ez);
      else if (dimension == TE)
        interpolateTimeDomainFieldCentralE_TE(E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, i, j, k, &Ex, &Ey,
                                              &Ez);
      else
        interpolateTimeDomainFieldCentralE_TM(E.xy, E.xz, E.yx, E.yz, E.zx, E.zy, i, j, k, &Ex, &Ey,
                                              &Ez);
      if (dimension == THREE)
        if (J_tot == 0) {
          interpolateTimeDomainFieldCentralH_2Dy(H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, i, j, k, &Hx,
                                                 &Hy, &Hz);
        } else if (intmethod == 1)
          interpolateTimeDomainFieldCentralH(H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, i, j, k, &Hx, &Hy,
                                             &Hz);
        else
          interpolateTimeDomainFieldCentralHBandLimited(H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, i, j, k,
                                                        &Hx, &Hy, &Hz);
      else if (dimension == TE)
        interpolateTimeDomainFieldCentralH_TE(H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, i, j, k, &Hx, &Hy,
                                              &Hz);
      else
        interpolateTimeDomainFieldCentralH_TM(H.xy, H.xz, H.yx, H.yz, H.zx, H.zy, i, j, k, &Hx, &Hy,
                                              &Hz);

      update_EH(EHr, EHi, vindex, campssample.components.index(FieldComponents::Ex), cphaseTermE, Ex);
      update_EH(EHr, EHi, vindex, campssample.components.index(FieldComponents::Hx), cphaseTermH, Hx);
      update_EH(EHr, EHi, vindex, campssample.components.index(FieldComponents::Ey), cphaseTermE, Ey);
      update_EH(EHr, EHi, vindex, campssample.components.index(FieldComponents::Hy), cphaseTermH, Hy);
      update_EH(EHr, EHi, vindex, campssample.components.index(FieldComponents::Ez), cphaseTermE, Ez);
      update_EH(EHr, EHi, vindex, campssample.components.index(FieldComponents::Hz), cphaseTermH, Hz);
    }
  }//end parallel region
}


void update_EH(double **EHr, double **EHi, int vindex, int idx, complex<double> &phase_term, double &value){

  if (idx >= 0) {
    auto tmp = value * phase_term; //exp(phaseTermE * I) * 1./((double) Nt);
    EHr[idx][vindex] += real(tmp);
    EHi[idx][vindex] += imag(tmp);
  }
}


void extractPhasorsSurfaceNoInterpolation(double **surface_EHr, double **surface_EHi,
                                          ElectricSplitField &E, MagneticSplitField &H,
                                          int **surface_vertices, int n_surface_vertices, int n,
                                          double omega, int Nt, int J_tot, SimulationParameters &params) {
  int vindex;
  double Ex, Ey, Ez, Hx, Hy, Hz;
  complex<double> phaseTermE, phaseTermH, subResultE, subResultH, cphaseTermE, cphaseTermH;

  phaseTermE = fmod(omega * ((double) n) * params.dt, 2 * dcpi);
  phaseTermH = fmod(omega * ((double) n + 0.5) * params.dt, 2 * dcpi);

  cphaseTermH = exp(phaseTermH * I) * 1. / ((double) Nt);
  cphaseTermE = exp(phaseTermE * I) * 1. / ((double) Nt);

  //loop over every vertex in the list
#pragma omp parallel default(shared) private(Ex, Ey, Ez, Hx, Hy, Hz, phaseTermE, phaseTermH,       \
                                             subResultE, subResultH, vindex)
  {
#pragma omp for
    for (vindex = 0; vindex < n_surface_vertices; vindex++) {
      //    fprintf(stderr,"vindex: %d: (%d %d %d)\n",vindex,surface_vertices[0][vindex],surface_vertices[1][vindex],surface_vertices[2][vindex]);
      Ex = E.xy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           E.xz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Ey = E.yx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           E.yz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Ez = E.zx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           E.zy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];

      Hx = H.xy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           H.xz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Hy = H.yx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           H.yz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Hz = H.zx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           H.zy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];

      /*Ex and Hx*/
      subResultH = Hx * cphaseTermH;//exp(phaseTermH * I) * 1./((double) Nt);
      subResultE = Ex * cphaseTermE;//exp(phaseTermE * I) * 1./((double) Nt);

      //now update the master array
      surface_EHr[0][vindex] = surface_EHr[0][vindex] + real(subResultE);
      surface_EHi[0][vindex] = surface_EHi[0][vindex] + imag(subResultE);

      surface_EHr[3][vindex] = surface_EHr[3][vindex] + real(subResultH);
      surface_EHi[3][vindex] = surface_EHi[3][vindex] + imag(subResultH);

      /*Ey and Hy*/
      subResultH = Hy * cphaseTermH;//exp(phaseTermH * I) * 1./((double) Nt);
      subResultE = Ey * cphaseTermE;//exp(phaseTermE * I) * 1./((double) Nt);

      //now update the master array
      surface_EHr[1][vindex] = surface_EHr[1][vindex] + real(subResultE);
      surface_EHi[1][vindex] = surface_EHi[1][vindex] + imag(subResultE);

      surface_EHr[4][vindex] = surface_EHr[4][vindex] + real(subResultH);
      surface_EHi[4][vindex] = surface_EHi[4][vindex] + imag(subResultH);


      /*Ez and Hz*/
      subResultH = Hz * cphaseTermH;//exp(phaseTermH * I) * 1./((double) Nt);
      subResultE = Ez * cphaseTermE;//exp(phaseTermE * I) * 1./((double) Nt);

      //now update the master array
      surface_EHr[2][vindex] = surface_EHr[2][vindex] + real(subResultE);
      surface_EHi[2][vindex] = surface_EHi[2][vindex] + imag(subResultE);

      surface_EHr[5][vindex] = surface_EHr[5][vindex] + real(subResultH);
      surface_EHi[5][vindex] = surface_EHi[5][vindex] + imag(subResultH);
    }
  }//end parallel region
}

void extractPhasorsPlane(double **iwave_lEx_Rbs, double **iwave_lEx_Ibs, double **iwave_lEy_Rbs,
                         double **iwave_lEy_Ibs, double **iwave_lHx_Rbs, double **iwave_lHx_Ibs,
                         double **iwave_lHy_Rbs, double **iwave_lHy_Ibs,
                         ElectricSplitField &E, MagneticSplitField &H, int I_tot, int J_tot, int K1, int n,
                         double omega, double dt, int Nt) {

  complex<double> phaseTerm = 0., subResult = 0.;


  phaseTerm = fmod(omega * ((double) n) * dt, 2 * dcpi);
  int i, j;

  for (j = 0; j < J_tot; j++)
    for (i = 0; i < (I_tot + 1); i++) {


      //Eyz
      subResult = (E.yz[K1][j][i] + E.yx[K1][j][i]) * exp(phaseTerm * I) * 1. / ((double) Nt);

      iwave_lEy_Rbs[j][i] = iwave_lEy_Rbs[j][i] + real(subResult);
      iwave_lEy_Ibs[j][i] = iwave_lEy_Ibs[j][i] + imag(subResult);

      //Hxz
      subResult = (H.xz[K1 - 1][j][i] + H.xy[K1][j][i]) * exp(phaseTerm * I) * 1. / ((double) Nt);

      iwave_lHx_Rbs[j][i] = iwave_lHx_Rbs[j][i] + real(subResult);
      iwave_lHx_Ibs[j][i] = iwave_lHx_Ibs[j][i] + imag(subResult);
    }

  for (j = 0; j < (J_tot + 1); j++)
    for (i = 0; i < I_tot; i++) {


      //Exz
      subResult = (E.xz[K1][j][i] + E.xy[K1][j][i]) * exp(phaseTerm * I) * 1. / ((double) Nt);

      iwave_lEx_Rbs[j][i] = iwave_lEx_Rbs[j][i] + real(subResult);
      iwave_lEx_Ibs[j][i] = iwave_lEx_Ibs[j][i] + imag(subResult);

      //Hyz
      subResult = (H.yz[K1 - 1][j][i] + H.yx[K1][j][i]) * exp(phaseTerm * I) * 1. / ((double) Nt);

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

double checkPhasorConvergence(ElectricField &E, ElectricField &E_copy) {

  double max_abs = 0., max_abs_diff = 0.;

  //find the largest maximum absolute value the largest difference (in absolute value) between phasors
  for (char c : {'x', 'y', 'z'})
    for (int k = 0; k < E.K_tot; k++)
      for (int j = 0; j < E.J_tot; j++)
        for (int i = 0; i < E.I_tot; i++){

            auto E_ijk = E.real[c][k][j][i] + I * E.imag[c][k][j][i];
            auto E_copy_ijk = E_copy.real[c][k][j][i] + I * E_copy.imag[c][k][j][i];

            max_abs = max(max_abs, abs(E_ijk));  // max(max_abs, |Re(E_x) + i Im(E_x)|)
            max_abs_diff = max(max_abs_diff, abs(E_ijk - E_copy_ijk));
          }

  return max_abs_diff / max_abs;
}

/*Copy the phasors from E to E_copy */
void copyPhasors(ElectricField &from, ElectricField &to, int nelements) {

  for (char c : {'x', 'y', 'z'}){
    memcpy(to.real[c], from.real[c], nelements * sizeof(double));
    memcpy(to.imag[c], from.imag[c], nelements * sizeof(double));
  }
}

/*Load up the output grid labels*/
void setGridLabels(GridLabels &input_labels, GridLabels &output_labels, int i_l, int i_u, int j_l,
                   int j_u, int k_l, int k_u) {
  //fprintf(stderr,"Entered: setGridLabels\n");
  //fprintf(stderr,"setGridLabels: %d,%d\n",j_u,j_l);
  for (int i = i_l; i <= i_u; i++) { output_labels.x[i - i_l] = input_labels.x[i]; }
  for (int j = j_l; j <= j_u; j++) { output_labels.y[j - j_l] = input_labels.y[j]; }
  for (int k = k_l; k <= k_u; k++) { output_labels.z[k - k_l] = input_labels.z[k]; }
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

/*work out if we have a conductive background*/
bool is_conductive(const XYZVectors &rho, int I_tot, int J_tot, int K_tot) {

  for (int i = 0; i < (I_tot + 1); i++)
    if (fabs(rho.x[i]) > 1e-15) return true;
  for (int j = 0; j < (J_tot + 1); j++)
    if (fabs(rho.y[j]) > 1e-15) return true;
  for (int k = 0; k < (K_tot + 1); k++)
    if (fabs(rho.z[k]) > 1e-15) return true;

  return false;
}

/*work out if we have a dispersive background*/
bool is_dispersive_ml(const DispersiveMultiLayer &ml, int K_tot) {
  for (int i = 0; i < K_tot; i++)
    if (fabs(ml.gamma[i]) > 1e-15) return true;
  return false;
}
