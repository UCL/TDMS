/**
 * @file test_hdf5_reader.cpp
 * @author William Graham
 * @brief Unit tests for HDF5Reader I/O functions that are non-specific to a
 * TDMS class or datatype
 */
#include "hdf5_io/hdf5_reader.h"

#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

using namespace std;
using tdms_tests::uint16s_to_string;
using tdms_unit_test_data::struct_testdata;

TEST_CASE("Read from a MATLAB struct") {
  HDF5Reader MATFile(struct_testdata);

  SECTION("Read numeric scalars") {
    /* Read scalar values from the MATLAB struct into the array. */
    // Initialise with values distinct from the expected values
    double one_half = 0., unity = 0.;
    bool logical_read = false;
    // Read values
    MATFile.read_field_from_struct("example_struct", "double_half", &one_half);
    MATFile.read_field_from_struct("example_struct", "double_no_decimal",
                                   &unity);
    MATFile.read_field_from_struct("example_struct", "boolean", &logical_read);
    // Validate read in data
    REQUIRE(one_half == Catch::Approx(0.5));
    REQUIRE(unity == Catch::Approx(1.));
    REQUIRE(int(unity) == 1);
    REQUIRE(logical_read);
  }

  SECTION("Read array data") {
    /* Read in the character array data */
    // string field is set to "tdms". Note that MATLAB saves this as uint16s, so
    // we need to convert manually...
    {
      uint16_t read_uints16[4];
      MATFile.read_field_from_struct("example_struct", "string", read_uints16);
      string tdms = uint16s_to_string(read_uints16, 4);
      REQUIRE(tdms == "tdms");
    }
    // The uint 3*4*5 uint matrix contains only 1s
    {
      vector<uint8_t> uint_matrix(3 * 4 * 5, 0);
      MATFile.read_field_from_struct("example_struct", "uint_345",
                                     uint_matrix.data());
      bool all_values_unity = true;
      for (uint8_t &value : uint_matrix) {
        if (value != 1) { all_values_unity = false; }
      }
      REQUIRE(all_values_unity);
    }
    // The double 2*2 matrix contains 0.25, 0.5, 0.75, 1.
    {
      double two_by_two[4];
      MATFile.read_field_from_struct("example_struct", "double_22", two_by_two);
      REQUIRE(two_by_two[0] == Catch::Approx(0.25));
      REQUIRE(two_by_two[1] == Catch::Approx(0.75));
      REQUIRE(two_by_two[2] == Catch::Approx(0.5));
      REQUIRE(two_by_two[3] == Catch::Approx(1.0));
    }
    // The complex matrix is the Pauli-y matrix [0, -i; i, 0]
  }
}
