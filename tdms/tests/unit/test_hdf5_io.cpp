/**
 * @file test_hdf5_io.cpp
 * @brief Tests of the HDF5 file I/O functionality.
 */
#include "hdf5_io.h"

// std
#include <ctime>
#include <filesystem>
#include <random>

// external
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

// tdms
#include "arrays.h"
#include "unit_test_utils.h"

using tdms_tests::create_tmp_dir;// unit_test_utils.h

TEST_CASE("Wrong datatype passed to ijk.") {
  auto tmp = create_tmp_dir();
  HDF5Writer r(tmp.string() + "/why.h5");
  // r.data_dump("x");
}

TEST_CASE("Test file I/O construction/destruction.") {

  // test-case wide setup - temporary directory
  auto tmp = create_tmp_dir();

  SECTION("Check file creation.") {
    HDF5Writer f(tmp.string() + "/test_file_constructor.h5");
    CHECK(f.is_ok());
  }

  SECTION("Check all reasonable file extensions are OK.") {
    for (auto extension : {".hdf5", ".h5", ".mat"}) {
      {
        HDF5Writer fw(tmp.string() + "/test_file" + extension);
        CHECK(fw.is_ok());

      }// Destructor called as we leave scope.

      HDF5Reader fr(tmp.string() + "/test_file" + extension);
      CHECK(fr.is_ok());
    }
  }

  SECTION("Check can't open nonexistent file.") {
    CHECK_THROWS(HDF5Reader(tmp.string() + "/this_file_doesnt_exist.h5"));
  }

  SECTION("Check can't read nonexistent data.") {
    {
      HDF5Writer fw(tmp.string() + "/this_file_does_exist_but_is_empty.h5");
      CHECK(fw.is_ok());

    }// Destructor called as we leave scope.

    double data[1];
    HDF5Reader fr(tmp.string() + "/this_file_does_exist_but_is_empty.h5");
    CHECK_THROWS(fr.read("nonexistantdata", data));
  }

  // Normal operation: we should be able to create a file and write to it, then
  // read from it.
  SECTION("Check write then read.") {
    {
      HDF5Writer fw(tmp.string() + "/test_file_wr.h5");
      hsize_t dimensions[1] = {1};
      double writeme = 1337.;
      fw.write("testdata", &writeme, 1, dimensions);
      SPDLOG_DEBUG("Written data");

      CHECK(fw.is_ok());
      fw.ls();

    }// Destructor called as we leave scope.

    double data[1];
    HDF5Reader fr(tmp.string() + "/test_file_wr.h5");
    fr.read("testdata", data);
    SPDLOG_DEBUG("Have read {}!", data[0]);
  }

  SECTION("Check write then (overwrite) then read.") {
    // Create the file and write some data.
    {
      HDF5Writer f1(tmp.string() + "/test_file_wor.h5");
      hsize_t dimensions[1] = {1};
      double writeme = 12345;
      f1.write("testdata", &writeme, 1, dimensions);
      SPDLOG_DEBUG("Written first data");

      CHECK(f1.is_ok());

    }// Destructor called as we leave scope.

    // Overwrite the file and add some different data.
    {
      HDF5Writer f2(tmp.string() + "/test_file_wor.h5");
      hsize_t dimensions[1] = {1};
      double writeme = 54321.;
      f2.write("testdata2", &writeme, 1, dimensions);
      SPDLOG_DEBUG("Written second data");

      CHECK(f2.is_ok());

    }// destructor called as we leave scope

    // Now open the file with a Reader. The first data should not be there (and
    // should throw an exception). The second data should be there.
    HDF5Reader f3(tmp.string() + "/test_file_wor.h5");

    CHECK(f3.get_datanames().size() == 1);

    double data[1];
    HDF5Reader fr(tmp.string() + "/test_file_wor.h5");
    CHECK_THROWS(f3.read("testdata", data));

    f3.read("testdata2", data);
    CHECK(data[0] == Catch::Approx(54321.));
  }

  // teardown - remove temporary directory and all files
  SPDLOG_DEBUG("Removing temporary directory.");
  std::filesystem::remove_all(tmp);
}

TEST_CASE("Test w/r TDMS objects") {
  // test-case wide setup - temporary directory
  auto tmp = create_tmp_dir();

  SECTION("5-element 1D array") {
    double to_write[5] = {1 / 137.0, 3.0, 2.71215, 3.14159, 916.0};
    double read_back[5];
    {
      hsize_t dimensions[1] = {5};
      HDF5Writer f1(tmp.string() + "/five_elements.h5");
      f1.write("five_elements", to_write, 1, dimensions);
    }
    {
      HDF5Reader f2(tmp.string() + "/five_elements.h5");
      f2.read("five_elements", read_back);
    }
    for (unsigned int i = 0; i < 5; i++) {
      CHECK(to_write[i] == Catch::Approx(read_back[i]));
    }
  }

  SECTION("5-by-6 2D array") {
    SPDLOG_INFO("5-by-6 2D array");
    Matrix<double> counting_matrix(5, 6);
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 6; j++) { counting_matrix[i][j] = 6. * i + j; }
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
        SPDLOG_INFO("Checking {} == {}", counting_matrix[i][j],
                    read_back[i][j]);
        CHECK(counting_matrix[i][j] == Catch::Approx(read_back[i][j]));
      }
    }
  }

  // teardown - remove temporary directory and all files
  SPDLOG_DEBUG("Removing temporary directory.");
  std::filesystem::remove_all(tmp);
}
