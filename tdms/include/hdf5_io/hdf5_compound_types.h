#pragma once

#include <complex>
#include <vector>

#include <H5Cpp.h>

#include "globals.h"

namespace hdf5_MATLAB_complex {
struct MATLAB_complex {
  double real;
  double imag;

  operator std::complex<double>() const {
    return real + imag * tdms_math_constants::IMAGINARY_UNIT;
  }

  double norm() const { return real * real + imag * imag; }
  double abs() const { return sqrt(norm()); }
};

H5::CompType to_hdf5_CompType();

std::vector<std::complex<double>>
to_complex_vector(const std::vector<MATLAB_complex> &matlab_complex);
}// namespace hdf5_MATLAB_complex
