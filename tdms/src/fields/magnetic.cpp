#include "field.h"
#include "interpolate_Hfield.h"

using namespace std;


double MagneticField::phase(int n, double omega, double dt){
  return omega * ((double) n + 0.5) * dt;  // 0.5 added because it's known half a time step after E
}

double MagneticSplitField::interpolate_x_to_centre(int i, int j, int k) {
  
}
double MagneticSplitField::interpolate_y_to_centre(int i, int j, int k) {}
double MagneticSplitField::interpolate_z_to_centre(int i, int j, int k) {}