#pragma once

#include "mat_io.h"

/**
 * @brief Defines a cuboid by specifying the first and last Yee cells in each
axial direction that form part of the cube.
 *
 * @details For example, { 0, 5, 2, 7, 1, 10 } corresponds to the cuboid that
contains all Yee cells indexed (i,j,k) where;
 * 0 <= i <= 5,
 * 2 <= j <= 7,
 * 1 <= k <= 10.
 *
 * TODO: Check inclusivity of the inequalities above.
 */
struct Cuboid {
  int array[6] = {0, 0, 0, 0, 0, 0};

  int operator[](int index) const { return array[index]; }
};
