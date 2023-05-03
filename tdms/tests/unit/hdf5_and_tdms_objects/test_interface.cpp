/**
 * @file test_hdf5_io.cpp
 * @brief Tests of the HDF5 file I/O functionality.
 */
#include "hdf5_io.h"

// std
// #include <cstdlib>
// #include <ctime>
// #include <filesystem>
// #include <random>
#include <string>

// external
// #include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

// tdms
#include "interface.h"
#include "unit_test_utils.h"

using namespace std;
using tdms_unit_test_data::tdms_object_data;

/**
 * @brief Check that HDF5 can read an InterfaceComponent from a HDF5 file.
 *
 * .mat files save InterfaceComponents as datasets of the "interface" group.
 *
 * interface.I0 = [1 0]; %I0 in the I=1 plane, no source condition
 * interface.I1 = [4 1]; %I1 in the I=4 plane, source condition applied
 * interface.J0 = [2 0]; %J0 in the J=2 plane, no source condition
 * interface.J1 = [5 0]; %J1 in the J=5 plane, no source condition
 * interface.K0 = [3 1]; %K0 in the K=3 plane, source condition applied
 * interface.K1 = [6 1]; %K1 in the K=6 plane, source condition applied
 */
TEST_CASE("HDF5: Read InterfaceComponent") {
  HDF5Reader MATFile(tdms_object_data);

  InterfaceComponent MATFile.read("interface")
}
