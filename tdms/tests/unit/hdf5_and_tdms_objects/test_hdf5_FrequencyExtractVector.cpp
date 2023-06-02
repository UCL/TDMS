#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays/vector_typedefs.h"
#include "globals.h"
#include "hdf5_io/hdf5_reader.h"
#include "unit_test_utils.h"

using Catch::Approx;
using tdms_math_constants::DCPI;
using tdms_unit_test_data::struct_testdata;
using tdms_unit_test_data::tdms_object_data;

TEST_CASE("HDF5Reader: FrequencyExtractVector") {
  FrequencyExtractVector fev;
  // This will make the default-assigned value 1.
  double omega_an = 2. * DCPI;

  SECTION("Empty dataset") {
    HDF5Reader input_file(struct_testdata);
    // Read from the empty array
    SECTION("Empty array") { input_file.read(fev, omega_an, "empty_array"); }
    // Read from the empty struct
    SECTION("Empty struct/group") {
      input_file.read(fev, omega_an, "empty_struct");
    }

    // Given the dataset that was pointed to was empty, fev should have been
    // assigned a single element equal to omega_an / 2 / DCPI.
    REQUIRE(fev.size() == 1);
    REQUIRE(fev[0] == Approx(omega_an / 2. / DCPI));
  }

  SECTION("Populated dataset") {
    // Since we only need to read in a vector of doubles in this case, we will
    // recycle the phasorsurface field which has the same format as f_ex_vec.
    HDF5Reader input_file(tdms_object_data);

    input_file.read(fev, omega_an, "phasorsurface");
    // The values we expected to read
    std::vector<double> expected_read = {1., 4., 2., 5., 3., 6.};

    // Validate the entries and length of the vector
    REQUIRE(fev == expected_read);
  }
}
