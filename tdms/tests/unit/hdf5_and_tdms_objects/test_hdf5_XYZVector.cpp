/**
 * @file test_hdf5_Matrix.cpp
 * @brief Tests of the HDF5 file I/O functionality when reading/writing Matrix
 * objects.
 */
#include "hdf5_io/hdf5_reader.h"

#include <string>

// external
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

// tdms
#include "arrays/xyz_vector.h"
#include "unit_test_utils.h"

using Catch::Approx;
using std::string;
using tdms_unit_test_data::hdf5_test_file;
using tdms_unit_test_data::struct_testdata;

TEST_CASE("HDF5: Read XYZVector") {
  XYZVector xyz;
  string group_name = "XYZVector";
  string prefix = "xyz_";

  SECTION("Read from .mat file") {
    HDF5Reader mat_file(struct_testdata);
    mat_file.read(group_name, prefix, xyz);
  }
  SECTION("Read from .hdf5 file") {
    HDF5Reader hdf5_file(hdf5_test_file);
    hdf5_file.read(group_name, prefix, xyz);
  }

  // Check data read in was OK
  bool all_entries_read_correctly = true;
  for (int i = 0; i < 3; i++) {
    all_entries_read_correctly = all_entries_read_correctly &&
                                 xyz['x'][i] == Approx((i + 1) / 10.) &&
                                 xyz['y'][i] == Approx((i + 4) / 10.) &&
                                 xyz['z'][i] == Approx((i + 7) / 10.);
  }
  REQUIRE(all_entries_read_correctly);
}
