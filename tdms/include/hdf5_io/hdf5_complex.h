#pragma once

#include <complex>

#include <H5Cpp.h>

#include "globals.h"

namespace hdf5_MATLAB_complex {
struct MATLAB_complex {
  double real;
  double imag;

  operator std::complex<double>() const {
    return real + imag * tdms_math_constants::IMAGINARY_UNIT;
  }
};

H5::CompType hdf5_CompType() {
  H5::CompType hdf5_MATLAB_complex(sizeof(MATLAB_complex));

  hdf5_MATLAB_complex.insertMember("real", HOFFSET(MATLAB_complex, real),
                                   H5::PredType::NATIVE_DOUBLE);
  hdf5_MATLAB_complex.insertMember("imag", HOFFSET(MATLAB_complex, imag),
                                   H5::PredType::NATIVE_DOUBLE);

  return hdf5_MATLAB_complex;
}
}// namespace hdf5_MATLAB_complex
