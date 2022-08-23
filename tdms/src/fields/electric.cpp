#include "globals.h"
#include "field.h"

using namespace std;


void ElectricField::add_to_angular_norm(int n, int Nt, SimulationParameters &params) {
  angular_norm += phasor_norm(ft, n, params.omega_an, params.dt, Nt);
}

double ElectricField::phase(int n, double omega, double dt){
  return omega * ((double) n + 1) * dt;
}
