#include "simulation_manager/simulation_manager.h"

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

void SimulationManager::extract_phasors(int &dft_counter, int tind) {
  // if we are not performing a complete run, these steps are skipped
  if (inputs.params.run_mode != RunMode::complete) { return; }

  // extract phasors if we are running a steady-state simulation
  if ((inputs.params.source_mode == SourceMode::steadystate) &&
      inputs.params.exphasorsvolume) {
    outputs.E.set_phasors(inputs.E_s, dft_counter - 1, inputs.params.omega_an,
                          inputs.params.dt, inputs.Nsteps);
    outputs.H.set_phasors(inputs.H_s, dft_counter, inputs.params.omega_an,
                          inputs.params.dt, inputs.Nsteps);
    // if we are additionally extracting surface phasors
    if (inputs.params.exphasorssurface) {
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        outputs.surface_phasors.extractPhasorsSurface(
                ifx, inputs.E_s, inputs.H_s, dft_counter,
                inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.Nsteps, inputs.params,
                inputs.params.intphasorssurface);
      }
      dft_counter++;
    }
  } else {
    // Running a pulsed simulation

    // Extract the phasors in the volume
    if (inputs.params.exphasorsvolume) {
      // compute volume phasors
      if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
      if ((tind - inputs.params.start_tind) % inputs.params.Np == 0) {
        outputs.E.set_phasors(inputs.E_s, tind - 1, inputs.params.omega_an,
                              inputs.params.dt, inputs.params.Npe);
        outputs.H.set_phasors(inputs.H_s, tind, inputs.params.omega_an,
                              inputs.params.dt, inputs.params.Npe);
      }
      if (TIME_EXEC) { timers.click_timer(TimersTrackingLoop::INTERNAL); }
    }
    // Extract phasors on the user-defined surface
    if ((inputs.params.exphasorssurface) &&
        ((tind - inputs.params.start_tind) % inputs.params.Np == 0)) {
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        outputs.surface_phasors.extractPhasorsSurface(
                ifx, inputs.E_s, inputs.H_s, tind,
                inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params.Npe,
                inputs.params, inputs.params.intphasorssurface);
      }
    }
    // Extract phasors at the user-defined vertices
    if (outputs.vertex_phasors.there_are_vertices_to_extract_at() &&
        ((tind - inputs.params.start_tind) % inputs.params.Np == 0)) {
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        outputs.vertex_phasors.extractPhasorsVertices(
                ifx, inputs.E_s, inputs.H_s, tind,
                inputs.f_ex_vec[ifx] * 2 * DCPI, inputs.params);
      }
    }
  }
}
