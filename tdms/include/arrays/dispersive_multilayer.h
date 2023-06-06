#pragma once

#include <vector>

#include "arrays/xyz_vector.h"

struct DispersiveMultiLayer {
public:
  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<double> gamma;
  XYZVector kappa;
  XYZVector sigma;

  /**
   * @brief Determines whether the (background) medium is dispersive
   *
   * @param near_zero_tolerance Tolerance for non-zero gamma (attenuation)
   * values
   * @return true Background is dispersive
   * @return false Background is not dispersive
   */
  bool is_dispersive(double near_zero_tolerance = 1e-15) const;
};
