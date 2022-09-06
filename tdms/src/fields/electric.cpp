#include "globals.h"
#include "field.h"

using namespace std;


double ElectricField::phase(int n, double omega, double dt){
  return omega * ((double) n + 1) * dt;
}
