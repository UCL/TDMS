/**
 * @file test_hdf5_io.cpp
 * @brief Tests of the HDF5 file I/O functionality.
 */
#include "hdf5_io/hdf5_reader.h"

// external
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
 * In the .mat file, the planes are saved as:
 * interface.I0 = [1 0]; %I0 in the I=1 plane, no source condition
 * interface.I1 = [4 1]; %I1 in the I=4 plane, source condition applied
 * interface.J0 = [2 0]; %J0 in the J=2 plane, no source condition
 * interface.J1 = [5 0]; %J1 in the J=5 plane, no source condition
 * interface.K0 = [3 1]; %K0 in the K=3 plane, source condition applied
 * interface.K1 = [6 1]; %K1 in the K=6 plane, source condition applied
 *
 * Do not forget the index offset when initialising from MATLAB indices! The
 * first element in each of these arrays should be offset by -1 upon being read
 * in.
 *
 * We will explicitly compare bools to false in what follows to make it
 * explicitly clear that we are testing that the values read into our object
 * _match_ those we expect from the data file.
 */
TEST_CASE("HDF5: Read InterfaceComponent") {
  tdms_unit_test_data::skip_if_missing(tdms_object_data);
  HDF5Reader MATFile(tdms_object_data);

  SECTION("Read into existing InterfaceComponent") {
    InterfaceComponent I0, J0, K0;

    MATFile.read("I0", &I0);
    MATFile.read("J0", &J0);
    MATFile.read("K0", &K0);

    bool I0_correct = (I0.index == 0) && (I0.apply == false);
    bool J0_correct = (J0.index == 1) && (J0.apply == false);
    bool K0_correct = (K0.index == 2) && (K0.apply == true);
    CHECK(I0_correct);
    CHECK(J0_correct);
    CHECK(K0_correct);
  }

  SECTION("Return InterfaceComponent object") {
    InterfaceComponent I1 = MATFile.read("I1");
    InterfaceComponent J1 = MATFile.read("J1");
    InterfaceComponent K1 = MATFile.read("K1");

    bool I1_correct = (I1.index == 3) && (I1.apply == true);
    bool J1_correct = (J1.index == 4) && (J1.apply == false);
    bool K1_correct = (K1.index == 5) && (K1.apply == true);
    CHECK(I1_correct);
    CHECK(J1_correct);
    CHECK(K1_correct);
  }
}
