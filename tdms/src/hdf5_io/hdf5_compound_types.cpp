#include "hdf5_io/hdf5_compound_types.h"

using namespace std;

H5::CompType hdf5_MATLAB_complex::to_hdf5_CompType() {
  H5::CompType hdf5_MATLAB_complex(sizeof(MATLAB_complex));

  hdf5_MATLAB_complex.insertMember("real", HOFFSET(MATLAB_complex, real),
                                   H5::PredType::NATIVE_DOUBLE);
  hdf5_MATLAB_complex.insertMember("imag", HOFFSET(MATLAB_complex, imag),
                                   H5::PredType::NATIVE_DOUBLE);

  return hdf5_MATLAB_complex;
}

vector<complex<double>> hdf5_MATLAB_complex::to_complex_vector(
        const vector<MATLAB_complex> &matlab_complex) {
  vector<complex<double>> std_complex(matlab_complex.size());

  // Use the defined casting method to assign to a vector of
  // std::complex<double>
  for (unsigned int i = 0; i < matlab_complex.size(); i++) {
    std_complex[i] = matlab_complex[i];
  }

  return std_complex;
}
