/**
 * @file test_hdf5_reader.cpp
 * @author William Graham
 * @brief Unit tests for HDF5Reader I/O functions that are non-specific to a
 * TDMS class or datatype
 */
#include "hdf5_io/hdf5_reader.h"

#include <stdint.h>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <spdlog/spdlog.h>

#include "unit_test_utils.h"

using namespace std;
using tdms_tests::uint16s_to_string;
using tdms_unit_test_data::struct_testdata, tdms_unit_test_data::hdf5_test_file;

TEST_CASE("HDF5Reader::read() [.mat files]") {
  HDF5Reader MATFile(struct_testdata);

  SECTION("Read numeric scalars") {
    /* Read scalar values from the MATLAB struct into the array. */
    // Initialise with values distinct from the expected values
    double one_half = 0., unity = 0.;
    bool logical_read = false;
    // Read values
    MATFile.read("example_struct/double_half", &one_half);
    MATFile.read("example_struct/double_no_decimal", &unity);
    MATFile.read("example_struct/boolean", &logical_read);
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
      MATFile.read("example_struct/string", read_uints16);
      string tdms = uint16s_to_string(read_uints16, 4);
      REQUIRE(tdms == "tdms");
    }
    // The uint 3*4*5 uint matrix contains only 1s
    {
      vector<uint8_t> uint_matrix(3 * 4 * 5, 0);
      MATFile.read("example_struct/uint_345", uint_matrix);
      bool all_values_unity = true;
      for (uint8_t &value : uint_matrix) {
        if (value != 1) { all_values_unity = false; }
      }
      REQUIRE(all_values_unity);
    }
    // The double 2*2 matrix contains 0.25, 0.5, 0.75, 1.
    {
      double two_by_two[4];
      MATFile.read("example_struct/double_22", two_by_two);
      REQUIRE(two_by_two[0] == Catch::Approx(0.25));
      REQUIRE(two_by_two[1] == Catch::Approx(0.75));
      REQUIRE(two_by_two[2] == Catch::Approx(0.5));
      REQUIRE(two_by_two[3] == Catch::Approx(1.0));
    }
  }
}

/** @brief Test the performance of read(), on both MATLAB files
 * and HDF5 files */
TEST_CASE("HDF5Reader::read() [std::vector]") {
  vector<double> read_buffer;
  read_buffer.reserve(12);
  // Used to check that data has been read in correctly
  bool entries_read_correctly = true;

  SECTION(".mat files") {
    HDF5Reader Hfile(struct_testdata);

    SECTION("Vector [int32]") {
      // We do an extra loop over the entries here to ensure that our final
      // check confirms that both the int data and the doubles that they were
      // cast to are correct
      vector<int> int_buffer(12);
      Hfile.read("read_in_test/vector", int_buffer);
      for (int i = 0; i < 12; i++) {
        entries_read_correctly = entries_read_correctly && int_buffer[i] == i;
        read_buffer[i] = (double) int_buffer[i];
      }
    }
    SECTION("Matrix [double]") {
      Hfile.read("read_in_test/matrix", read_buffer);
    }
    SECTION("Tensor [double]") {
      Hfile.read("read_in_test/tensor", read_buffer);
    }
  }

  SECTION(".hdf5 files") {
    HDF5Reader Hfile(hdf5_test_file);

    // h5py saves int dtype at 64-bit integers
    SECTION("Vector [int64]") {
      // We do an extra loop over the entries here to ensure that our final
      // check confirms that both the int data and the doubles that they were
      // cast to are correct
      vector<int64_t> int_buffer(12);
      Hfile.read("read_in_test/vector_int", int_buffer);
      for (int i = 0; i < 12; i++) {
        entries_read_correctly = entries_read_correctly && int_buffer[i] == i;
        read_buffer[i] = (double) int_buffer[i];
      }
    }
    SECTION("Matrix [double]") {
      Hfile.read("read_in_test/matrix_double", read_buffer);
    }
    SECTION("Tensor [double]") {
      Hfile.read("read_in_test/tensor_double", read_buffer);
    }
  }

  for (int i = 0; i < 12; i++) {
    entries_read_correctly =
            entries_read_correctly && read_buffer[i] == Catch::Approx(i);
  }
  REQUIRE(entries_read_correctly);
}
