#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "hdf5_io/hdf5_reader.h"
#include "unit_test_utils.h"

using tdms_unit_test_data::hdf5_test_file;
using tdms_unit_test_data::struct_testdata;

// HDF5Base cannot be instantiated, so we use HDF5Reader just to avoid
// accidental writes.

TEST_CASE("HDF5Base::path_exists()") {
  HDF5Reader matlab_file(struct_testdata);
  HDF5Reader hdf5_file(hdf5_test_file);

  SECTION("Paths that exist") {
    spdlog::info("HDF5Base::path_exists() [Valid paths]");
    // Point to a dataset that exists
    REQUIRE(matlab_file.path_exists("empty_array"));
    // Point to a group that exists
    REQUIRE(matlab_file.path_exists("empty_struct"));
  }

  SECTION("Enforce object types") {
    spdlog::info("HDF5Base::path_exists() [Forced types]");
    // read_in_test/vector is a dataset
    REQUIRE(matlab_file.path_exists("read_in_test/vector", H5I_DATASET));
    REQUIRE(!matlab_file.path_exists("read_in_test/vector", H5I_GROUP));

    // read_in_test is a group
    REQUIRE(hdf5_file.path_exists("read_in_test", H5I_GROUP));
    REQUIRE(!hdf5_file.path_exists("read_in_test", H5I_ATTR));
  }

  SECTION("Request object types") {
    spdlog::info("HDF5Base::path_exists() [Request types]");
    H5I_type_t dataset;
    H5I_type_t group;

    // Point to a dataset
    REQUIRE(matlab_file.path_exists("read_in_test/vector", &dataset));
    REQUIRE(dataset == H5I_DATASET);

    // Point to a group
    REQUIRE(hdf5_file.path_exists("read_in_test", &group));
    REQUIRE(group == H5I_GROUP);
  }

  SECTION("Paths that do not exist") {
    spdlog::info("HDF5Base::path_exists() [Invalid paths]");
    REQUIRE(!hdf5_file.path_exists("this_path_is_invalid"));
    REQUIRE(!matlab_file.path_exists("read_in_test/this_dataset_doesnt_exist"));
  }
}

TEST_CASE("HDF5Base: flagged_MATLAB_empty()") {
  HDF5Reader matlab_file(struct_testdata);
  HDF5Reader hdf5_file(hdf5_test_file);

  SECTION("Point to empty objects") {
    spdlog::info("HDF5Base::flagged_MATLAB_empty() [Empty objects]");
    // Point to empty dataset (array)
    REQUIRE(matlab_file.flagged_MATLAB_empty("empty_array"));
    // Point to empty group (struct)
    REQUIRE(matlab_file.flagged_MATLAB_empty("empty_struct"));
  }

  SECTION("Point to non-empty objects") {
    spdlog::info("HDF5Base::flagged_MATLAB_empty() [Non-empty objects]");
    // Point to nonempty dataset (array)
    REQUIRE(!matlab_file.flagged_MATLAB_empty("nonempty_dataset"));
    // Point to nonempty group (struct)
    REQUIRE(!matlab_file.flagged_MATLAB_empty("example_struct"));
    // Point to a nonempty dataset within a group
    REQUIRE(!matlab_file.flagged_MATLAB_empty("read_in_test/vector"));
    // Point to nonempty group not created by MATLAB (so the attribute should
    // not be flagged)
    REQUIRE(!hdf5_file.flagged_MATLAB_empty("read_in_test"));
  }

  SECTION("Exception cases") {
    spdlog::info("HDF5Base::flagged_MATLAB_empty() [Exceptions]");
    // Object not present in the file
    REQUIRE_THROWS_AS(
            matlab_file.flagged_MATLAB_empty("no_object_with_this_name"),
            std::runtime_error);
    REQUIRE_THROWS_AS(
            hdf5_file.flagged_MATLAB_empty("no_object_with_this_name"),
            std::runtime_error);

    // Object is not a group or dataset
    REQUIRE_THROWS_AS(hdf5_file.flagged_MATLAB_empty("file_attribute"),
                      std::runtime_error);
  }

  SECTION("False positives") {
    spdlog::info("HDF5Base::flagged_MATLAB_empty() [False positives]");
    // Object has an attribute MATLAB_empty, but still contains data
    REQUIRE(hdf5_file.flagged_MATLAB_empty("flag_as_empty"));
    // Object is an empty array, has the MATLAB_empty flag, but this is set to
    // false and thus it should come back as not being empty
    REQUIRE(!hdf5_file.flagged_MATLAB_empty("not_marked_empty"));
  }
}
