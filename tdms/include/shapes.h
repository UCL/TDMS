#pragma once

#include "mat_io.h"

class Cuboid {
private:
  /* The indices of the first and last Yee cells in each axial direction that
  are encompassed by the cuboid. For example, { 0, 5, 2, 7, 1, 10 } corresponds
  to the cuboid that contains all Yee cells indexed (i,j,k) where 0 <= i <= 5, 2
  <= j <= 7, 1 <= k <= 10.

  TODO: Check inclusivity of the inequalities above.
  */
  int array[6] = {0, 0, 0, 0, 0, 6};

public:
  void initialise(const mxArray *ptr, int J_tot);

  inline int operator[](int value) const { return array[value]; };
};
