#include <stdexcept>
#include "dimensions.h"


using namespace std;

int Dimensions::operator[](int value) const {
  switch (value) {
    case 0: return i;
    case 1: return j;
    case 2: return k;
    default: throw runtime_error("Have no element " + to_string(value));
  }
}

Dimensions::Dimensions(const mxArray *ptr) {

  auto n_dims = mxGetNumberOfDimensions(ptr);
  auto raw_dims = mxGetDimensions(ptr);

  if (n_dims > 3){
    throw runtime_error("Cannot initialise more than 3D");
  }

  switch (n_dims) {
    case 3: k = raw_dims[2];
    case 2: j = raw_dims[1];
    case 1: i = raw_dims[0];
    default: break;
  }
}
