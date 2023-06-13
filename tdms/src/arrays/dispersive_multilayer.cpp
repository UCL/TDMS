#include "arrays/dispersive_multilayer.h"

bool DispersiveMultiLayer::is_dispersive(double near_zero_tolerance) const {
  for (double gamma_val : gamma) {
    if (fabs(gamma_val) > near_zero_tolerance) {
      // non-zero attenuation constant of a Yee cell implies media is dispersive
      return true;
    }
  }
  return false;
}
