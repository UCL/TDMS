/**
 * @file dispersive_multilayer.h
 * @author William Graham
 * @brief Defines DispersiveMultilayer, a structure that holds dispersion
 * constants for the material being simulated.
 */
#pragma once

#include <vector>

#include "arrays/xyz_vector.h"

/**
 * @brief Stores dispersive constants that describe the properties of the
 *        multilayered structure being simulated.
 */
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
