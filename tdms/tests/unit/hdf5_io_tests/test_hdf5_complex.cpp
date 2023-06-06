#include <complex>
#include <vector>

#include <H5Cpp.h>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "hdf5_io/hdf5_complex.h"
#include "hdf5_io/hdf5_reader.h"
#include "unit_test_utils.h"

using namespace std;
using hdf5_MATLAB_complex::MATLAB_complex;
using tdms_unit_test_data::struct_testdata;

TEST_CASE("Compound data do gooder?") {
  spdlog::info("Declare things");
  vector<MATLAB_complex> data_vector;
  HDF5Reader matlab_file(struct_testdata);

  spdlog::info("Do a read?");
  matlab_file.read("example_struct/complex_22", data_vector,
                   hdf5_MATLAB_complex::hdf5_CompType());

  for (const MATLAB_complex &c : data_vector) {
    spdlog::info("This member is (real) {} + (imag) {}i", c.real, c.imag);
    complex<double> std_c = c;
    spdlog::info("Cast has returned: {} + {}i", std_c.real(), std_c.imag());
  }
}
