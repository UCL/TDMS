/**
 * @file main.cpp
 * @brief The main function. Launches TDMS.
 */
#include <vector>
#include <string>

#include <spdlog/spdlog.h>

#include "openandorder.h"
#include "input_output_names.h"
#include "input_matrices.h"
#include "iterator.h"
#include "mat_io.h"

using namespace tdms_matrix_names;

int main(int nargs, char *argv[]){

  // Set the logging level with a compile-time #define for debugging
  #if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_DEBUG
    spdlog::set_level(spdlog::level::debug);
  #elif SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_INFO
    spdlog::set_level(spdlog::level::info);
  #endif

  mxArray *plhs[NOUTMATRICES_PASSED];
  InputMatrices matrix_inputs;

  auto args = ArgumentParser::parse_args(nargs, argv);
  check_files_can_be_accessed(args);

  //now it is safe to use matlab routines to open the file and order the matrices
  if (!args.has_grid_filename()) {
    matrix_inputs.set_from_input_file(args.input_filename());
  } else {
    matrix_inputs.set_from_input_file(args.input_filename(), args.grid_filename());
  }

  // decide which derivative method to use (PSTD or FDTD)
  SolverMethod solver_method = PseudoSpectral; // default
  if (args.finite_difference())
    solver_method = SolverMethod::FiniteDifference;

  // decide whether to toggle off the band-limited interpolation methods
  PreferredInterpolationMethods preferred_interpolation_methods =
          PreferredInterpolationMethods::BandLimited;// default
  if (args.cubic_interpolation()) {
    preferred_interpolation_methods = PreferredInterpolationMethods::Cubic;
  }

  // now run the time propagation code
  spdlog::info("\n==== Executing simulation ==== \n");
  execute_simulation(NOUTMATRICES_PASSED, (mxArray **) plhs, NMATRICES,
                     matrix_inputs, solver_method, preferred_interpolation_methods);
  spdlog::info("\n==== Recieved output data from the simulation ====\n");

  if (!args.have_flag("-m")) {//prints vertices and facets
    saveoutput(plhs, matricestosave_all, outputmatrices_all, NOUTMATRICES_WRITE_ALL,
               args.output_filename());
  } else {// minimise the file size by not printing vertices and facets
    saveoutput(plhs, matricestosave, outputmatrices, NOUTMATRICES_WRITE,
               args.output_filename());
  }

  return 0;
}
