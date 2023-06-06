#pragma once

#include <vector>

#include "globals.h"

struct XYZVector {
  std::vector<double> x = {};
  std::vector<double> y = {};
  std::vector<double> z = {};
  //   /**
  //    * Set the pointer for one of the vectors in this collection with a name
  //    of
  //    * c
  //    * @param c Character labeling the vector
  //    * @param ptr Pointer to assign
  //    */
  //   void set_ptr(char c, double *ptr);
  //   /**
  //    * Set the pointer for one of the vectors in this collection with a name
  //    of
  //    * c
  //    * @param d AxialDirection labeling the vector
  //    * @param ptr Pointer to assign
  //    */
  //   void set_ptr(AxialDirection d, double *ptr);

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
