#include "iterator.h"

#include <omp.h>
#include <spdlog/spdlog.h>

#include "iterator_class.h"
#include "matrix.h"

using namespace std;

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
  // SETUP SIMULATION AND REPORT CONFIG OPTIONS
  spdlog::info("== Setting up simulation ==\n");

  // report number of threads that will be used
  spdlog::info("Using {0:d} OMP threads", omp_get_max_threads());

  // validate that we recieved the correct number of inputs and outputs to this function
  if (nrhs != 49) { throw runtime_error("Expected 49 inputs. Had " + to_string(nrhs)); }
  if (nlhs != 31) { throw runtime_error("31 outputs required. Had " + to_string(nlhs)); }

  // read in the information from the input files and command-line, and setup the variables to be used in the main loop
  Iterator main_loop(in_matrices, solver_method, preferred_interpolation_methods);

  // link loop variables to the output pointers in plhs
  main_loop.link_fields_and_labels(plhs);
  main_loop.link_id(plhs);
  main_loop.link_fdtd_phasor_arrays(plhs);
  main_loop.link_fieldsample(plhs);
  main_loop.link_vertex_phasors(plhs);

  // Perpare the parameters used in the phasor convergence procedure
  main_loop.prepare_phasor_convergence_proceedure();

  // Perform the j-loop optimisation, if possible
  main_loop.optimise_loops_if_possible();

  // [MAIN LOOP] RUN TDMS SIMULATION

  spdlog::debug("== Starting main loop ==\n");
  main_loop.run_main_loop();

  // POST-LOOP PROCESSING

  // normalise the phasors in the volume (if we are extracting them)
  main_loop.normalise_field_volumes();
  // normalise the phasors on the surface wrt the {E,H}-phasor-norms
  main_loop.normalise_surface_phasors();
  // normalise the phasors at the user-requested vertices
  main_loop.normalise_vertex_phasors();
  // normalise the Id output array data
  main_loop.normalise_Id_arrays();

  // OUTPUT ASSIGNMENT

  // find the maximum absolute value of residual field in the grid
  double maxfield = main_loop.compute_max_split_field_value();
  // place this in the output
  int ndims = 2;
  int dims[2] = {1, 1};
  plhs[25] = mxCreateNumericArray(ndims, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  *mxGetPr((mxArray *) plhs[25]) = maxfield;

  // write the interpolated fields and the corresponding gridlabels
  if (main_loop.params.run_mode == RunMode::complete && main_loop.params.exphasorsvolume) {
    // do some writing - STILL NEED TO CLEANUP
    main_loop.initialise_output_labels_from_input_labels();

    // now interpolate over the extracted phasors
    main_loop.interpolate_field_values(plhs);

    //now set up the output grid labels for the interpolated fields
    int Ex_label_dims[2] = {1, main_loop.E.I_tot - 3};
    plhs[19] = mxCreateNumericArray(2, (const mwSize *) Ex_label_dims, mxDOUBLE_CLASS, mxREAL);

    int Ey_label_dims[2] = {1, max(1, main_loop.E.J_tot - 3)};
    plhs[20] = mxCreateNumericArray(2, (const mwSize *) Ey_label_dims, mxDOUBLE_CLASS, mxREAL);

    int Ez_label_dims[2] = {1, 1};
    if (main_loop.params.dimension == Dimension::THREE) { Ez_label_dims[1] = main_loop.E.K_tot - 3; }
    plhs[21] = mxCreateNumericArray(2, (const mwSize *) Ez_label_dims, mxDOUBLE_CLASS, mxREAL);

    // write the interpolated coordinates to the output
    main_loop.write_interpolated_gridlabels(plhs);
  } else {
    // we do not want to write the interpolated fields, set output to be empty arrays
    int emptydims[2] = {0, 0};
    for (int emptyloop = 13; emptyloop <= 21; emptyloop++) {
      plhs[emptyloop] =
              mxCreateNumericArray(2, (const mwSize *) emptydims, mxDOUBLE_CLASS, mxCOMPLEX);
    }
  }

  // Now export 3 matrices: a vertex list, a matrix of complex amplitudes at these vertices, and a list of facets
  if (main_loop.params.exphasorssurface && main_loop.params.run_mode == RunMode::complete) {
    main_loop.regenerate_mesh_for_facets(plhs);
  } else {
    // set outputs as empty arrays, user did not request this information
    int emptydims[2] = {0, 0};
    plhs[22] = mxCreateNumericArray(2, (const mwSize *) emptydims, mxDOUBLE_CLASS, mxREAL);
    plhs[23] = mxCreateNumericArray(2, (const mwSize *) emptydims, mxDOUBLE_CLASS, mxREAL);
    plhs[24] = mxCreateNumericArray(2, (const mwSize *) emptydims, mxDOUBLE_CLASS, mxREAL);
  }

  // END OUTPUT PROCESSING AND WRITING

  // ~Iterator now handles tear down and cleanup. plhs is returned to main(). plhs pointers are preserved, intermediate pointers to the data location of the saved arrays should be cleaned up by ~Iterator. Hanging MATLAB memory should be free'd by the appropriate classes going out of scope.
}
