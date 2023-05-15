#include "simulation_manager/simulation_manager.h"

void SimulationManager::end_of_iteration_steps(double &time_of_last_log,
                                               unsigned int tind,
                                               ElectricField &E_copy) {
  // If enough time has passed since the last write to the log
  if ((((double) time(NULL)) - time_of_last_log) > 1) {
    // Compute and print max residual field
    double maxfield = std::max(inputs.E_s.largest_field_value(),
                               inputs.H_s.largest_field_value());
    spdlog::info("Iterating: tind = {0:d}, maxfield = {1:e}", tind, maxfield);
    // Update time of last log to now
    time_of_last_log = double(time(NULL));
  }

  // Report on possible phasor convergence failure, if:
  if ((inputs.params.source_mode ==
       SourceMode::steadystate) &&// running a steady state simulation
      (tind ==
       (inputs.params.Nt -
        1)) &&// have performed the allowed number of algorithm iterations
      (inputs.params.run_mode ==
       RunMode::complete) &&        // complete run mode requested
      inputs.params.exphasorsvolume)// extracting phasors in the volume
  {
    spdlog::info("Iteration limit reached (no convergence): setting output "
                 "fields to last "
                 "complete DFT");
    outputs.E.set_values_from(E_copy);
  }

  // Export the field if we are at a suitable iteration number
  if (inputs.params.has_tdfdir && (tind % inputs.params.Np) == 0) {
    spdlog::info("Exporting field...");
    inputs.ex_td_field_exporter.export_field(inputs.E_s, inputs.skip_tdf, tind);
  }
}
