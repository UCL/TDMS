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
#include "output_matrices.h"

using namespace tdms_matrix_names;

int main(int nargs, char *argv[]){

  // Set the logging level with a compile-time #define for debugging
  #if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_DEBUG
    spdlog::set_level(spdlog::level::debug);
  #elif SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_INFO
    spdlog::set_level(spdlog::level::info);
  #endif

  InputMatrices matrix_inputs;
  OutputMatrices outputs;

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
  execute_simulation(NOUTMATRICES_PASSED, outputs, NMATRICES,
                     matrix_inputs, solver_method, preferred_interpolation_methods);

  // save the outputs, possibly in compressed format
  outputs.save_outputs(args.output_filename(), args.have_flag("-m"));

  return 0;
}
