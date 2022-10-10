#include <algorithm>
#include "shapes.h"


using namespace std;

void Cuboid::initialise(const mxArray *ptr, int J_tot) {

  auto ndims = mxGetNumberOfDimensions(ptr);
  auto dims = mxGetDimensions((mxArray *) ptr);

  if (ndims != 2) {
    throw runtime_error("expected phasorsurface to be a vector of length 6");
  }
  if (dims[0] != 1 || dims[1] != 6) {
    throw runtime_error("expected phasorsurface to be a vector of length 6");
  }

  for (int i = 0; i < 6; i++) {
    array[i] = max((int) *(mxGetPr(ptr) + i) - 1, // Change from MATLAB -> C indexing
                   0);
  }
  if (J_tot == 0 && array[2] != array[3]) {
      throw runtime_error("When doing a 2D simulation, J0 should equal J1 in phasorsurface.");
  }
}
