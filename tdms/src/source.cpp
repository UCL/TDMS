#include "source.h"

#include <iostream>
#include <spdlog/spdlog.h>
#include <stdexcept>

#include "dimensions.h"
#include "matlabio.h"


using namespace std;

Source::Source(const mxArray *ptr, int dim1, int dim2,
               const std::string &name) {

  if (mxIsEmpty(ptr)) {
    spdlog::info("{} is empty", name);
  } else {
    // fetch dimensions of the input array, and check they are what we expect
    auto dims = Dimensions(ptr);

    if (dims.are_1d()) {
      throw runtime_error(name + " should be 3- or 2-dimensional");
    }
    if (dims.are_2d()) { dim2 = 0; }
    if (!(dims[0] == 8 && dims[1] == dim1 && dims[2] == dim2)) {
      cerr << name << " has incorrect size" << endl;
    }

    // check that complex data has been passed
    if (!mxIsComplex(ptr)) {
      throw runtime_error(
              name +
              " should be complex, use a call of "
              "complex(real(Isource),imag(Isource)) in matlab if necessary");
    }

    // cast MATLAB arrays to C++ datatypes
    real = cast_matlab_3D_array(mxGetPr(ptr), dims[0], dims[1], dims[2]);
    imag = cast_matlab_3D_array(mxGetPi(ptr), dims[0], dims[1], dims[2]);
    // flag Source as non-empty
    no_data_stored = false;
  }
}
