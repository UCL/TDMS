/**
 * @file test_hdf5_Matrix.cpp
 * @brief Tests of the HDF5 file I/O functionality when reading/writing Matrix
 * objects.
 */
#include "hdf5_io/hdf5_reader.h"
#include "hdf5_io/hdf5_writer.h"

#include <filesystem>

// external
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

// tdms
#include "arrays.h"
#include "unit_test_utils.h"

using std::filesystem::remove_all;
using tdms_tests::create_tmp_dir;

TEST_CASE("HDF5: Read/Write Matrix") {
  // test-case wide setup - temporary directory
  auto tmp = create_tmp_dir();

  SECTION("5-by-6 2D array [double]") {
    SPDLOG_INFO("5-by-6 2D array");
    Matrix<double> counting_matrix(5, 6);
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 6; j++) { counting_matrix(i, j) = 6. * i + j; }
    }
    Matrix<double> read_back;

    {
      HDF5Writer f1(tmp.string() + "/five-by-six.h5");
      f1.write("five-by-six", counting_matrix);
    }
    {
      SPDLOG_DEBUG("About to read...");
      HDF5Reader f2(tmp.string() + "/five-by-six.h5");
      f2.read("five-by-six", read_back);
    }

    for (unsigned int i = 0; i < 5; i++) {
      for (unsigned int j = 0; j < 6; j++) {
        SPDLOG_INFO("Checking {} == {}", counting_matrix(i, j),
                    read_back(i, j));
        CHECK(counting_matrix(i, j) == Catch::Approx(read_back(i, j)));
      }
    }
  }

  // teardown - remove temporary directory and all files
  SPDLOG_DEBUG("Removing temporary directory.");
  remove_all(tmp);
}
