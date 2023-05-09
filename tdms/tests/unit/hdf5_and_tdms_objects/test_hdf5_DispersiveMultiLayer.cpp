#include "hdf5_io/hdf5_reader.h"

#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

using tdms_unit_test_data::tdms_object_data;

TEST_CASE("HDF5: Read DispersiveMultiLayer") {
  HDF5Reader MATFile(tdms_object_data);
  // read from dispersive_aux group
  DispersiveMultiLayer dml;

  MATFile.read(&dml);
}
