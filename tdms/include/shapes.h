#include "matlabio.h"


class Cuboid{
private:
  int array[6] = {0, 0, 0, 0, 0, 6};

public:
  void initialise(const mxArray *ptr, int J_tot);

  inline int operator[] (int value) const { return array[value]; };
};
