/**
 * @file main.cpp
 * @brief The main function. Launches TDMS.
 */
#include "argument_parser.h"
#include "input_file_settings.h"
#include "input_matrices.h"
#include "simulation_manager/simulation_manager.h"

using namespace tdms_matrix_names;

int main(int nargs, char *argv[]) {

// Set the logging level with a compile-time #define for debugging
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_DEBUG
  spdlog::set_level(spdlog::level::debug);
#elif SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_INFO
  spdlog::set_level(spdlog::level::info);
#endif

  InputMatrices matrix_inputs;

  auto args = ArgumentParser::parse_args(nargs, argv);
  args.check_files_can_be_accessed();

  // now it is safe to use matlab routines to open the file and order the
  // matrices
  if (!args.has_grid_filename()) {
    matrix_inputs.set_from_input_file(args.input_filename());
  } else {
    matrix_inputs.set_from_input_file(args.input_filename(),
                                      args.grid_filename());
  }

  InputFileSettings settings(args.input_filename());
  fprintf(stdout, "usecd: %d, intmethod %d\n", settings.solver(),
          settings.interpolation());

  // decide which derivative method to use (PSTD or FDTD)
  SolverMethod solver_method = SolverMethod::PseudoSpectral;// default
  if (args.finite_difference()) solver_method = SolverMethod::FiniteDifference;

  // decide whether to toggle off the band-limited interpolation methods
  InterpolationMethod preferred_interpolation_methods =
          InterpolationMethod::BandLimited;// default
  if (args.cubic_interpolation()) {
    preferred_interpolation_methods = InterpolationMethod::Cubic;
  }

  // Handles the running of the simulation, given the inputs to the executable.
  SimulationManager simulation(matrix_inputs, solver_method,
                               preferred_interpolation_methods);

  // now run the time propagation code
  simulation.execute();

  // perform post-loop processes on the outputs
  simulation.post_loop_processing();

  // save the outputs, possibly in compressed format
  simulation.write_outputs_to_file(args.output_filename(),
                                   args.have_flag("-m"));

  return 0;
}
