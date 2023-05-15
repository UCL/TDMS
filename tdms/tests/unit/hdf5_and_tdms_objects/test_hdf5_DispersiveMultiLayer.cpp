#include "hdf5_io/hdf5_reader.h"

#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

using Catch::Approx;
using namespace std;
using tdms_unit_test_data::cant_find_test_data_message,
        tdms_unit_test_data::tdms_object_data;

TEST_CASE("HDF5: Read DispersiveMultiLayer") {
  if (!std::filesystem::exists(tdms_object_data))
    SKIP(cant_find_test_data_message);
  HDF5Reader MATFile(tdms_object_data);
  // read from dispersive_aux group
  DispersiveMultiLayer dml;

  SECTION("Correct data") {
    vector<double> consecutive_integers(10);
    for (int i = 0; i < 10; i++) { consecutive_integers[i] = (double) i; }

    MATFile.read(&dml);

    // Assert correct data - each entry should just be the integers 0->9
    // inclusive
    bool correct_data_read = true;
    for (int i = 0; i < 10; i++) {
      correct_data_read = correct_data_read && dml.alpha[i] == Approx(i) &&
                          dml.beta[i] == Approx(i) && dml.gamma[i] == Approx(i);
    }
    correct_data_read = correct_data_read &&
                        equal(dml.kappa.x.begin(), dml.kappa.x.end(),
                              consecutive_integers.begin()) &&
                        equal(dml.kappa.y.begin(), dml.kappa.y.end(),
                              consecutive_integers.begin()) &&
                        equal(dml.kappa.z.begin(), dml.kappa.z.end(),
                              consecutive_integers.begin()) &&
                        equal(dml.sigma.x.begin(), dml.sigma.x.end(),
                              consecutive_integers.begin()) &&
                        equal(dml.sigma.y.begin(), dml.sigma.y.end(),
                              consecutive_integers.begin()) &&
                        equal(dml.sigma.z.begin(), dml.sigma.z.end(),
                              consecutive_integers.begin());

    REQUIRE(correct_data_read);
  }
}
