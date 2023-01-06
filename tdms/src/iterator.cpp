#include "iterator.h"

#include <omp.h>
#include <spdlog/spdlog.h>

#include "iterator.h"
#include "matrix.h"

using namespace std;

void execute_simulation(int nlhs, mxArray *plhs[], int nrhs, InputMatrices in_matrices,
                        SolverMethod solver_method,
                        PreferredInterpolationMethods preferred_interpolation_methods) {
  // SETUP SIMULATION AND REPORT CONFIG OPTIONS
  spdlog::info("\n== Setting up simulation ==\n");

  // report number of threads that will be used
  spdlog::info("Using {0:d} OMP threads", omp_get_max_threads());

  // validate that we recieved the correct number of inputs and outputs to this function
  if (nrhs != 49) { throw runtime_error("Expected 49 inputs. Had " + to_string(nrhs)); }
  if (nlhs != 31) { throw runtime_error("31 outputs required. Had " + to_string(nlhs)); }

  // read in the information from the input files and command-line, and setup the variables to be used in the main loop
  spdlog::info("Setting up loop variables");
  Iterator main_loop(in_matrices, solver_method, preferred_interpolation_methods);

  // link loop variables to the output pointers in plhs
  spdlog::info("Linking to output pointers [pre-main-loop]");
  main_loop.link_fields_and_labels(plhs);
  main_loop.link_id(plhs);
  main_loop.link_fdtd_phasor_arrays(plhs);
  main_loop.link_fieldsample(plhs);
  main_loop.link_vertex_phasors(plhs);

  // Perpare the parameters used in the phasor convergence procedure
  spdlog::info("Preparing phasor convergence proceedure");
  main_loop.prepare_phasor_convergence_proceedure();

  // Perform the j-loop optimisation, if possible
  spdlog::info("Optimising j-loop");
  main_loop.optimise_loops_if_possible();

  // [MAIN LOOP] RUN TDMS SIMULATION
  spdlog::info("\n== Starting main loop ==\n");

  main_loop.run_main_loop();

  // POST-LOOP PROCESSING
  spdlog::info("\n== Processing outputs ==\n");

  // normalise the phasors in the volume (if we are extracting them)
  main_loop.normalise_field_volumes();
  // normalise the phasors on the surface wrt the {E,H}-phasor-norms
  main_loop.normalise_surface_phasors();
  // normalise the phasors at the user-requested vertices
  main_loop.normalise_vertex_phasors();
  // normalise the Id output array data
  main_loop.normalise_Id_arrays();

  // OUTPUT ASSIGNMENT
  spdlog::info("\n== Assigning outputs ==\n");

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
  spdlog::info("End simulation, beginning tear-down");

  // ~Iterator now handles tear down and cleanup. plhs is returned to main(). plhs pointers are preserved, intermediate pointers to the data location of the saved arrays should be cleaned up by ~Iterator. Hanging MATLAB memory should be free'd by the appropriate classes going out of scope.
}
