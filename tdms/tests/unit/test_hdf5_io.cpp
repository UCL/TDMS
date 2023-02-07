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
#include "unit_test_utils.h"

using tdms_tests::create_tmp_dir;// unit_test_utils.h

TEST_CASE("Test file I/O construction/destruction.") {
  // setup - temporary directory
  auto tmp = create_tmp_dir();

  SECTION("Check file creation.") {
    HDF5Writer f(tmp.string() + "/test_file_constructor.h5");
    CHECK(f.is_ok());
  }// destructor called as we leave scope

  SECTION("Check all reasonable file extensions are OK.") {
    for (auto extension : {".hdf5", ".h5", ".mat"}) {
      {
        HDF5Writer fw(tmp.string() + "/test_file" + extension);
        CHECK(fw.is_ok());
      }// destructor called as we leave scope

      HDF5Reader fr(tmp.string() + "/test_file" + extension);
      CHECK(fr.is_ok());
    }
  }

  SECTION("Check can't open nonexistent file.") {
    CHECK_THROWS(HDF5Reader(tmp.string() + "/this_file_doesnt_exist.h5"));
  }

  // Normal operation: we should be able to create a file and write to it, then read from it.
  SECTION("Check write then read.") {
    // create a file
    {
      HDF5Writer fw(tmp.string() + "/test_file_wr.h5");
      hsize_t dimensions[1] = {1};
      double writeme = 1337.;
      fw.write("testdata", &writeme, 1, dimensions);
      spdlog::debug("Written data");

      CHECK(fw.is_ok());
      fw.ls();

    }// destructor called as we leave scope

    double data[1];
    HDF5Reader fr(tmp.string() + "/test_file_wr.h5");
    fr.read("testdata", data);
    spdlog::debug("Have read {}!", data[0]);
  }

  /*
  SECTION("Check write then (overwrite) then read.") {

    // Create the file and write some data.
    {
      HDF5Writer f1(tmp.string() + "/test_file_wwr.h5");
      //f1.write();
      CHECK(f1.is_ok());

    }// destructor called as we leave scope

    // Overwrite the file and add some different data.
    {
      HDF5Writer f2(tmp.string() + "/test_file_wwr.h5");
      //f2.write();
      CHECK(f2.is_ok());

    }// destructor called as we leave scope

    // Now open the file in readonly. The first data should not be there (and
    // should throw an exception). The second data should be there.
    HDF5Reader f3(tmp.string() + "/test_file_wwr.h5");
    //f3.read();

    //CHECK_NOTHROW(f3.read());
    //CHECK_THROWS(f3.read());
  }
  */

  // teardown - remove temporary directory and all files
  SPDLOG_DEBUG("Removing temporary directory.");
  std::filesystem::remove_all(tmp);
}
