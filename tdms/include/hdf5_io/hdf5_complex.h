#pragma once

#include <vector>

#include <H5Cpp.h>

typedef struct MATLAB_complex {
  std::vector<double> real;
  hvl_t dreal_handle;
  std::vector<double> imag;
  hvl_t dimag_handle;
};

H5::CompType hdf5_MATLAB_complex(sizeof(MATLAB_complex));
auto double_type = H5::PredType::NATIVE_DOUBLE;
auto var_double_type = H5::VarLenType(&double_type);

hdf5_MATLAB_complex.insertMember("real", HOFFSET(MATLAB_complex, dreal_handle),
                                 var_double_type);
hdf5_MATLAB_complex.insertMember("imag", HOFFSET(MATLAB_complex, dimag_handle),
                                 var_double_type);
