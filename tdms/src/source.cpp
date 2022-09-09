#include <stdexcept>
#include <iostream>
#include "source.h"
#include "dimensions.h"


using namespace std;

Source::Source(const mxArray *ptr, int dim1, int dim2, const std::string &name){

  if (mxIsEmpty(ptr)) {
    cerr << name << " is empty" << endl;
  } else {
    auto dims = Dimensions(ptr);

    if (dims.are_1d()){
      throw runtime_error(name+" should be 3- or 2-dimensional");
    }
    if (dims.are_2d()){
      dim2 = 0;
    }
    if (!(dims[0] == 8 && dims[1] == dim1 && dims[2] == dim2)){
      cerr << name << " has incorrect size" << endl;
    }
    if (!mxIsComplex(ptr)){
      throw runtime_error(name+" should be complex, use a call of "
                          "complex(real(Isource),imag(Isource)) in matlab if necessary");
    }
    real = cast_matlab_3D_array(mxGetPr(ptr), dims[0], dims[1], dims[2]);
    imag = cast_matlab_3D_array(mxGetPi(ptr), dims[0], dims[1], dims[2]);
  }
}
