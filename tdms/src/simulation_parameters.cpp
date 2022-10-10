#include <iostream>
#include <stdexcept>
#include "simulation_parameters.h"


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
    dimension = Dimension::TE;
  } else {
    dimension = Dimension::TM;
  }
}

void SimulationParameters::set_phasorinc(const double* vector) {
  phasorinc.x = (int) vector[0];
  phasorinc.y = (int) vector[1];
  phasorinc.z = (int) vector[2];
}

void SimulationParameters::set_Np(FrequencyExtractVector &f_ex_vec) {

  double f_max = f_ex_vec.max();
  Np = (int) floor(1. / (2.5 * dt * f_max));

  //calculate Npe, the temporal DFT will be evaluated whenever tind increments by Npe
  for (unsigned int tind = start_tind; tind < Nt; tind++){
    if ((tind - start_tind) % Np == 0) Npe++;
  }
  cerr << "Np=" << Np << " Nt=" << Nt << " Npe=" << Npe
       <<  " f_max=" << f_max << " Npraw=" << 2.5 * dt * f_max << endl;
}
