#include "hdf5_io/hdf5_reader.h"

// external
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

// tdms
#include "unit_test_utils.h"

using tdms_unit_test_data::tdms_object_data;

TEST_CASE("HDF5: Read FrequencyVector") {
  HDF5Reader MATFile(tdms_object_data);

  // FrequencyVectors are read from the structure array "f_vec", and consist of
  // real-valued arrays of the same length "fx_vec" and "fy_vec"
  SECTION("Read into existing FrequencyVectors object") {
    FrequencyVectors f_vec;
    MATFile.read(&f_vec);

    // confirm values are assigned correctly
  }

  SECTION("Return FrequencyVectors object") {
    FrequencyVectors f_vec = MATFile.read();

    // confirm values are assigned correctly
  }
}
