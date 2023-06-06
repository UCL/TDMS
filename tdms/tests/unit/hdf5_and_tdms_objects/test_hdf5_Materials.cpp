#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays/material_collections.h"
#include "hdf5_io/hdf5_reader.h"
#include "unit_test_utils.h"

using Catch::Approx;
using std::runtime_error;
using tdms_unit_test_data::hdf5_test_file;
using tdms_unit_test_data::tdms_object_data;

TEST_CASE("CMaterial") {
  spdlog::info("CMaterial");
  CMaterial c_material;
  const std::vector<double> expected_read = {0., 1., 2., 3., 4.};

  SECTION(".mat input") {
    HDF5Reader file(tdms_object_data);

    // Check error is thrown if we point to a group we can't read from
    REQUIRE_THROWS_AS(file.read(c_material, "Dmaterial"), runtime_error);

    file.read(c_material);
  }
  SECTION(".hdf5 input") {
    HDF5Reader file(hdf5_test_file);
    file.read(c_material);
  }

  bool data_read_correctly = true;
  for (const char &xyz : {'x', 'y', 'z'}) {
    data_read_correctly = data_read_correctly &&
                          c_material.a[xyz] == expected_read &&
                          c_material.b[xyz] == expected_read &&
                          c_material.c[xyz] == expected_read;
  }
  REQUIRE(data_read_correctly);
}

TEST_CASE("CCollection") {
  spdlog::info("CCollection");
  CCollection c_collection;

  // This is the size of the vectors we expect to read in if data was present
  const int expected_nonempty_size = 5;
  // If we only read in a group with 6 attributes, the .c member should be
  // empty. We use this variable to record whether this is the case or not.
  int expected_c_member_size;

  SECTION("Read with 9 elements") {
    HDF5Reader file(tdms_object_data);
    file.read(c_collection, "Cmaterial");

    REQUIRE(c_collection.is_disp_ml == true);
    // In this case, the data is setup so that the .c member has the same sizes
    // as the other .a and .b members
    expected_c_member_size = expected_nonempty_size;
  }
  SECTION("Read with 6 elements") {
    HDF5Reader file(hdf5_test_file);
    file.read(c_collection);

    REQUIRE(c_collection.is_disp_ml == false);
    // In this case, the .c member should remain empty
    expected_c_member_size = 0;
  }

  bool correct_sizes = true;
  for (const char &xyz : {'x', 'y', 'z'}) {
    correct_sizes = correct_sizes &&
                    c_collection.a[xyz].size() == expected_nonempty_size &&
                    c_collection.b[xyz].size() == expected_nonempty_size &&
                    c_collection.c[xyz].size() == expected_c_member_size;
  }
  REQUIRE(correct_sizes);
}

TEST_CASE("DMaterial and DCollection") {
  spdlog::info("DMaterial and DCollection");
  DMaterial d_material;
  DCollection d_collection;

  // This is the data we expect to be read into each member
  const std::vector<double> expected_read = {5., 6., 7., 8., 9.};

  SECTION(".mat") {
    HDF5Reader file(tdms_object_data);

    // We can't read from a group that doesn't contain exactly 6 members
    REQUIRE_THROWS_AS(file.read(d_material, "Cmaterial"), runtime_error);

    file.read(d_material);
    // DCollections should be able to read from Dmaterials as they're the same
    // structure
    file.read(d_collection, "Dmaterial");
  }

  SECTION(".hdf5") {
    HDF5Reader file(hdf5_test_file);

    file.read(d_material);
    // DCollections should be able to read from DMaterials
    // as they're the same structure
    file.read(d_collection, "Dmaterial");
  }

  bool correctly_read = true;
  for (const char &xyz : {'x', 'y', 'z'}) {
    correctly_read = correctly_read && d_material.a[xyz] == expected_read &&
                     d_material.b[xyz] == d_material.a[xyz] &&
                     d_material.a[xyz] == d_collection.a[xyz] &&
                     d_material.b[xyz] == d_collection.b[xyz];
  }
  REQUIRE(correctly_read);
}
