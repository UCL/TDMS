#include "hdf5_io/hdf5_reader.h"

// external
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

// tdms
#include "unit_test_utils.h"

using Catch::Approx;
using tdms_unit_test_data::tdms_object_data;

inline int EXPECTED_VEC_SIZE = 4;

/**
 * @brief Test the reading of frequency vectors from the input file, and failure
 * cases when they are of the incorrect dimension.
 *
 * f_vec is a group whose data consists of:
 * fx_vec: double[4] ROW vector whose entries are 1/4, 2/4, 3/4, 4/4
 * fy_vec: double[4] COL vector whose entries are -1/4, -2/4, -3/4, -4/4
 *
 * f_vec_bad is a group whose data consists of:
 * fx_vec: double[2,2] of zero entries
 * fy_vec: double[1,2] of zero entries.
 */
TEST_CASE("HDF5: Read FrequencyVector") {
  HDF5Reader MATFile(tdms_object_data);

  // FrequencyVectors are automatically read from the structure array "f_vec",
  // and consist of real-valued arrays of the same length "fx_vec" and "fy_vec"
  FrequencyVectors f_vec;

  SECTION("Read into existing FrequencyVectors object") {
    MATFile.read(&f_vec);
  }

  SECTION("Return FrequencyVectors object") { f_vec = MATFile.read(); }

  // Confirm expected sizes
  CHECK(f_vec.x.size() == EXPECTED_VEC_SIZE);
  CHECK(f_vec.y.size() == EXPECTED_VEC_SIZE);
  // Confirm correct assignment of entries
  bool x_entries_correct = true;
  bool y_entries_correct = true;
  for (int i = 0; i < EXPECTED_VEC_SIZE; i++) {
    x_entries_correct =
            x_entries_correct && (f_vec.x[i] == Approx((i + 1) / 4.));
    y_entries_correct =
            y_entries_correct && (f_vec.y[i] == Approx(-(i + 1) / 4.));
  }
  REQUIRE(x_entries_correct);
  REQUIRE(y_entries_correct);
}
