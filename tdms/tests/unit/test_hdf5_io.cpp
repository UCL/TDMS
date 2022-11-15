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

using tdms_tests::create_tmp_dir; // unit_test_utils.h

TEST_CASE("Test file I/O construction/destruction.") {
  // setup - temporary directory
  auto tmp = create_tmp_dir();

  SECTION("Check file I/O modes") {
    // check possible file i/o modes: create file, then open same file readonly, then open
    // rw, then create it again (even though it already exists)
    for (auto i : {HDF5FileMode::OVERWRITE, HDF5FileMode::READONLY, HDF5FileMode::READWRITE,
                   HDF5FileMode::OVERWRITE}) {
      SPDLOG_DEBUG("Constructing in mode {}", i);
      HDF5File fi(tmp.string() + "/test_file_constructor.h5", static_cast<HDF5FileMode>(i));
      CHECK(fi.isOK());
    }// destructor called as we leave scope
  }

  SECTION("Check all reasonable file extensions are OK") {
    for (auto extension : {".hdf5", ".h5", ".mat"}) {
      HDF5File fi(tmp.string() + "/test_file" + extension, OVERWRITE);
      CHECK(fi.isOK());
    }
  }

  SECTION("Check can't open nonexistent file") {
    CHECK_THROWS(HDF5File(tmp.string() + "/this_file_doesnt_exist.h5", HDF5FileMode::READONLY));
  }

  // We should get an exception if we try to write to a readonly file.
  SECTION("Check can't open file readonly and write to it") {
    // create a file
    {
      HDF5File f1(tmp.string() + "/test_file_read_only.h5", HDF5FileMode::OVERWRITE);
    }// destructor called as we leave scope

    HDF5File f2(tmp.string() + "/test_file_read_only.h5", HDF5FileMode::READONLY);
    CHECK_THROWS(f2.write());
  }

  // Normal operation: we should be able to create a file and write to it, then read from it.
  SECTION("Check write then read.") {
    // create a file
    {
      HDF5File f1(tmp.string() + "/test_file_wr.h5", HDF5FileMode::OVERWRITE);
      f1.write();
      CHECK(f1.isOK());

    }// destructor called as we leave scope

    HDF5File f2(tmp.string() + "/test_file_wr.h5", HDF5FileMode::READONLY);
    f2.read();
  }

  SECTION("Check write then (overwrite) then read.") {

    // Create the file and write some data.
    {
      HDF5File f1(tmp.string() + "/test_file_wwr.h5", HDF5FileMode::OVERWRITE);
      f1.write();
      CHECK(f1.isOK());

    }// destructor called as we leave scope

    // Overwrite the file and add some different data.
    {
      HDF5File f2(tmp.string() + "/test_file_wwr.h5", HDF5FileMode::OVERWRITE);
      f2.write();
      CHECK(f2.isOK());

    }// destructor called as we leave scope

    // Now open the file in readonly. The first data should not be there (and
    // should throw an exception). The second data should be there.
    HDF5File f3(tmp.string() + "/test_file_wwr.h5", HDF5FileMode::READONLY);
    f3.read();

    CHECK_NOTHROW(f3.read());
    CHECK_THROWS(f3.read());
  }

  // teardown - remove temporary directory and all files
  SPDLOG_DEBUG("Removing temporary directory.");
  std::filesystem::remove_all(tmp);
}

TEST_CASE("Test call example sandbox function") {
  // setup - temporary directory
  auto tmp = create_tmp_dir();

  SPDLOG_INFO("Calling the example");
  example_hdf5();
  SPDLOG_INFO("Example has run");

  // teardown - remove temporary directory and all files
  SPDLOG_DEBUG("Removing temporary directory.");
  std::filesystem::remove_all(tmp);
}
