#include <catch2/catch_test_macros.hpp>

#include "hdf5_io/hdf5_reader.h"
#include "unit_test_utils.h"

using tdms_unit_test_data::hdf5_test_file;
using tdms_unit_test_data::struct_testdata;

// HDF5Base cannot be instantiated, so we use HDF5Reader just to avoid
// accidental writes.
TEST_CASE("HDF5Base: is_empty()") {
  HDF5Reader matlab_file(struct_testdata);
  HDF5Reader hdf5_file(hdf5_test_file);

  SECTION("Point to empty objects") {
    // Point to empty dataset (array)
    REQUIRE(matlab_file.has_MATLAB_empty("empty_array"));
    // Point to empty group (struct)
    REQUIRE(matlab_file.has_MATLAB_empty("empty_struct"));
  }

  SECTION("Point to non-empty objects") {
    // Point to nonempty dataset (array)
    REQUIRE(!matlab_file.has_MATLAB_empty("nonempty_dataset"));
    // Point to nonempty group (struct)
    REQUIRE(!matlab_file.has_MATLAB_empty("example_struct"));
    // Point to nonempty dataset not created by MATLAB (so the attribute should
    // not be flagged)
    REQUIRE(!hdf5_file.has_MATLAB_empty("read_in_test"));
  }

  SECTION("Exception cases") {
    // Object not present in the file
    REQUIRE_THROWS_AS(matlab_file.has_MATLAB_empty("no_object_with_this_name"),
                      std::runtime_error);
    REQUIRE_THROWS_AS(hdf5_file.has_MATLAB_empty("no_object_with_this_name"),
                      std::runtime_error);

    // Object is not a group or dataset
    REQUIRE_THROWS_AS(hdf5_file.has_MATLAB_empty("file_attribute"),
                      std::runtime_error);
  }
}
