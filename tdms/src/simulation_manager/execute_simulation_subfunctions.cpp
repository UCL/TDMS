#include "simulation_manager.h"

bool SimulationManager::check_phasor_convergence(int &dft_counter,
                                                 ElectricField &E_copy) {
  if ((dft_counter == inputs.Nsteps) &&
      (inputs.params.run_mode == RunMode::complete) &&
      (inputs.params.source_mode == SourceMode::steadystate) &&
      inputs.params.exphasorsvolume) {

    dft_counter = 0;

    double tol = outputs.E.normalised_difference(E_copy);
    if (tol < TOL) { return true; }// required accuracy obtained

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
  return false;
}
