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
