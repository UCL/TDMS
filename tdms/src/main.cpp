/**
 * @file main.cpp
 * @brief The main function. Launches TDMS.
 */
#include "argument_parser.h"
#include "simulation_manager/simulation_manager.h"

using namespace tdms_matrix_names;

int main(int nargs, char *argv[]) {

// Set the logging level with a compile-time #define for debugging
#if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_DEBUG
  spdlog::set_level(spdlog::level::debug);
#elif SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_INFO
  spdlog::set_level(spdlog::level::info);
#endif

  auto args = ArgumentParser::parse_args(nargs, argv);
  args.check_files_can_be_accessed();

  // Handles the running of the simulation, given the inputs to the executable.
  SimulationManager simulation =
          args.has_grid_filename() ? SimulationManager(args.input_filename(),
                                                       args.grid_filename())
                                   : SimulationManager(args.input_filename());

  // now run the time propagation code
  simulation.execute();

  // perform post-loop processes on the outputs
  simulation.post_loop_processing();

  // save the outputs, possibly in compressed format
  simulation.write_outputs_to_file(args.output_filename(),
                                   args.have_flag("-m"));

  return 0;
}
