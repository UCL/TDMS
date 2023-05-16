/**
 * @file test_hdf5_dimension.cpp
 * @author William Graham
 * @brief Unit tests for the H5Dimension object
 */
#include "hdf5_io/hdf5_dimension.h"
#include "hdf5_io/hdf5_reader.h"

#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

using namespace std;
using tdms_unit_test_data::struct_testdata;

TEST_CASE("Fetch dimensions correctly") {
  tdms_unit_test_data::skip_if_missing(struct_testdata);
  HDF5Reader MATFile(struct_testdata);

  SECTION("2D array") {
    H5Dimension two_by_two = MATFile.shape_of("example_struct", "double_22");
    bool two_by_two_dimensions_read_in =
            (two_by_two == vector<hsize_t>{2, 2}) && (!two_by_two.is_1D());
    REQUIRE(two_by_two_dimensions_read_in);
  }

  SECTION("3D array") {
    H5Dimension three_four_five =
            MATFile.shape_of("example_struct", "uint_345");
    bool three_four_five_read_in =
            (three_four_five == vector<hsize_t>{3, 4, 5}) &&
            (!three_four_five.is_1D());
  }

  SECTION("1D char array") {
    H5Dimension tdms_dims = MATFile.shape_of("example_struct", "string");
    bool tdms_dims_correct = (tdms_dims.is_1D()) && (tdms_dims.max_dim() == 4);
    REQUIRE(tdms_dims_correct);
  }
}
