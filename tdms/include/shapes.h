/**
 * @file shapes.h
 * @brief Classes describing geometric shapes.
 */
#pragma once

#include "mat_io.h"

/**
 * @brief Cuboid shape.
 * @details Used to hold the surface over which to extract phasors.
 */
class Cuboid {
private:
  int array[6] = {0, 0, 0, 0, 0, 6};//!< The surface.

public:
  /**
   * @brief Read in phasorsurface from a MATLAB array.
   *
   * @param ptr [in] A pointer to the MATLAB array.
   * @param J_tot [in] The total number of Yee cells in the y direction.
   */
  void initialise(const mxArray *ptr, int J_tot);

  /**
   * @brief Access the underlying data array.
   *
   * @param value The array index.
   * @return int The array value.
   */
  inline int operator[](int value) const { return array[value]; };
};
