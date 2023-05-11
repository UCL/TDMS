#include "simulation_parameters.h"

#include <iostream>
#include <stdexcept>

#include <spdlog/spdlog.h>

#include "arrays/incident_field.h"

using namespace std;

SimulationParameters::SimulationParameters() = default;

void SimulationParameters::set_source_mode(string mode_string) {

  auto s = move(mode_string);
  if (s == "steadystate") {
    source_mode = SourceMode::steadystate;
  } else if (s == "pulsed") {
    source_mode = SourceMode::pulsed;
  } else {
    throw runtime_error("value of sourcemode (" + s + ") is invalid");
  }
}

void SimulationParameters::set_run_mode(string mode_string) {

  auto s = move(mode_string);
  if (s == "complete") {
    run_mode = RunMode::complete;
  } else if (s == "analyse") {
    run_mode = RunMode::analyse;
  } else {
    throw runtime_error("value of runmode (" + s + ") is invalid");
  }
}

void SimulationParameters::set_dimension(string mode_string) {

  auto s = move(mode_string);
  if (s == "3") {
    dimension = Dimension::THREE;
  } else if (s == "TE") {
    dimension = Dimension::TRANSVERSE_ELECTRIC;
  } else {
    dimension = Dimension::TRANSVERSE_MAGNETIC;
  }
}

void SimulationParameters::set_spacing_stride(const double *vector) {
  spacing_stride.x = (int) vector[0];
  spacing_stride.y = (int) vector[1];
  spacing_stride.z = (int) vector[2];
}

void SimulationParameters::set_Np(FrequencyExtractVector &f_ex_vec) {

  double f_max = f_ex_vec.max();
  Np = (int) floor(1. / (2.5 * dt * f_max));

  // calculate Npe, the temporal DFT will be evaluated whenever tind increments
  // by Npe
  for (unsigned int tind = start_tind; tind < Nt; tind++) {
    if ((tind - start_tind) % Np == 0) Npe++;
  }
  spdlog::info("Np = {}, Nt = {}, Npe = {}, f_max = {}, Npraw = {}", Np, Nt,
               Npe, f_max, 2.5 * dt * f_max);
}

void SimulationParameters::unpack_from_input_matrices(
        InputMatrices in_matrices) {
  // determine if we have a dispersive medium or multilayer
  CCollection C(in_matrices["C"]);
  is_disp_ml = C.is_disp_ml;
  is_multilayer = C.is_multilayer;

  // set delta params
  delta.dx = *mxGetPr(ptr_to_vector_in(in_matrices["delta"], "x", "delta"));
  delta.dy = *mxGetPr(ptr_to_vector_in(in_matrices["delta"], "y", "delta"));
  delta.dz = *mxGetPr(ptr_to_vector_in(in_matrices["delta"], "z", "delta"));

  // unpack constants for the simulation
  omega_an = double_in(in_matrices["omega_an"], "omega_an");
  to_l = double_in(in_matrices["to_l"], "to_l");
  hwhm = double_in(in_matrices["hwhm"], "hwhm");
  pml.Dxl = int_cast_from_double_in(in_matrices["Dxl"], "Dxl");
  pml.Dxu = int_cast_from_double_in(in_matrices["Dxu"], "Dxu");
  pml.Dyl = int_cast_from_double_in(in_matrices["Dyl"], "Dyl");
  pml.Dyu = int_cast_from_double_in(in_matrices["Dyu"], "Dyu");
  pml.Dzl = int_cast_from_double_in(in_matrices["Dzl"], "Dzl");
  pml.Dzu = int_cast_from_double_in(in_matrices["Dzu"], "Dzu");
  Nt = int_cast_from_double_in(in_matrices["Nt"], "Nt");
  dt = double_in(in_matrices["dt"], "dt");
  start_tind = int_cast_from_double_in(in_matrices["tind"], "tind");

  // set the source and run mode
  set_source_mode(string_in(in_matrices["sourcemode"], "sourcemode"));
  set_run_mode(string_in(in_matrices["runmode"], "runmode"));

  // determine additional behaviour of this call to tdms
  exphasorsvolume = bool_cast_from_double_in(in_matrices["exphasorsvolume"],
                                             "exphasorsvolume");
  exphasorssurface = bool_cast_from_double_in(in_matrices["exphasorssurface"],
                                              "exphasorssurface");
  intphasorssurface = bool_cast_from_double_in(in_matrices["intphasorssurface"],
                                               "intphasorssurface");
  interp_mat_props =
          bool_cast_from_double_in(in_matrices["intmatprops"], "intmatprops");

  // set the stride and dimension of the simulation
  set_spacing_stride(mxGetPr(in_matrices["phasorinc"]));
  set_dimension(string_in(in_matrices["dimension"], "dimension"));

  // get exdetintegral if it is present
  if (!mxIsEmpty(in_matrices["exdetintegral"])) {
    exdetintegral = bool_cast_from_double_in(in_matrices["exdetintegral"],
                                             "exdetintegral");
  }

  // get air interface
  if (!mxIsEmpty(in_matrices["air_interface"])) {
    air_interface_present = true;
    air_interface = double_in(in_matrices["air_interface"], "air_interface");
  }

  // get the time-domain field
  IncidentField Ei(in_matrices["tdfield"]);
  exi_present = Ei.x.has_elements();
  eyi_present = Ei.y.has_elements();

  // set_Np -> do we actually USE this? f_ex_vec goes out of scope immediately
  // after...
  FrequencyExtractVector f_ex_vec(in_matrices["f_ex_vec"], omega_an);
  set_Np(f_ex_vec);
}
