#include "hdf5_io/hdf5_reader.h"

// external
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include <stdexcept>

#include "unit_test_utils.h"

using tdms_unit_test_data::tdms_object_data,
        tdms_unit_test_data::tdms_bad_object_data;

/**
 * @brief Test that a Cuboid object can be successfully read in.
 *
 * The (correct) Cuboid we have prepared is just the array {1, 4, 2, 5, 3, 6}.
 * The bad data provides us with a column vector of 7 elements.
 */
TEST_CASE("HDF5: Read Cuboid") {
  Cuboid cube;

  SECTION("Read into existing object") {
    HDF5Reader MATFile(tdms_object_data);
    MATFile.read(&cube);
    // Check expected values, noting the -1 offset that is applied because of
    // MATLAB indexing
    bool expected_values_recieved = cube[0] == 0 && cube[1] == 3 &&
                                    cube[2] == 1 && cube[3] == 4 &&
                                    cube[4] == 2 && cube[5] == 5;
    REQUIRE(expected_values_recieved);
  }

  SECTION("Throw error if too many elements provided") {
    HDF5Reader MATFile(tdms_bad_object_data);
    // Error should be thrown due to incorrect dimensions
    REQUIRE_THROWS_AS(MATFile.read(&cube), std::runtime_error);
  }
}
