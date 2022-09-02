#include "matlabio.h"
#include "field.h"
#include "globals.h"

using namespace std;

void Field::add_to_angular_norm(int n, int Nt, SimulationParameters &params) {
  angular_norm += phasor_norm(ft, n, params.omega_an, params.dt, Nt);
}
