/**
 * @file xyz_vector.h
 * @author William Graham
 * @brief Declares the XYZVector structure, which stores three vector members x,
 * y, z.
 */
#pragma once

#include <vector>

#include "globals.h"

/**
 * @brief Structure containing 3 std::vector<double> members: x, y, z.
 * @details Intended for use in instances where a quantity varies across the
 * three axial directions, but the variation in each direction is independent of
 * the other axial directions.
 *
 * Members can be referenced through the [] operator and passing the
 * corresponding single-letter character.
 */
struct XYZVector {
  std::vector<double> x = {};
  std::vector<double> y = {};
  std::vector<double> z = {};

  std::vector<double> &operator[](const char &direction) {
    switch (direction) {
      case 'x':
        return x;
      case 'y':
        return y;
      case 'z':
        return z;
      default:
        throw std::runtime_error("Component not recognised");
        break;
    }
  }

  /**
   * @brief Determines whether all elements in the x, y, or z vector are less
   * than a given value.
   *
   * @param comparison_value Value to compare elements to
   * @param component Vector to compare elements against; x, y, or z
   * @return true All elements are less than the comparison_value
   * @return false At least one element is not less than the comparison_value
   */
  bool all_elements_less_than(double comparison_value,
                              AxialDirection component) const;
  /**
   * @brief Determines whether all elements in the x, y, AND z vectors are
   * less than a given value.
   *
   * @param comparison_value Value to compare elements to
   * @return true All elements are less than the comparison_value
   * @return false At least one element is not less than the comparison_value
   */
  bool all_elements_less_than(double comparison_value) const;
};
