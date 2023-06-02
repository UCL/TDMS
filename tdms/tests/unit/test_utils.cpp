/**
 * @file test_utils.cpp
 * @brief Tests of the utility functions
 */
#include "utils.h"

#include <fstream>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

using Catch::Approx;
using std::vector;

TEST_CASE("Test nonexistant file") {
  REQUIRE_THROWS_AS(assert_can_open_file("DOES_NOT_EXIST.mat", "r"),
                    std::runtime_error);
}

TEST_CASE("Test existing file") {

  // create a oneline text file
  auto tmp = tdms_tests::create_tmp_dir();
  auto path = tmp / "DOES_EXIST.txt";
  std::ofstream dummyfile(path.c_str());
  dummyfile << "Hello" << std::endl;
  dummyfile.close();

  // should be able to open it
  REQUIRE_NOTHROW(assert_can_open_file(path.c_str(), "r"));
  REQUIRE_NOTHROW(assert_can_open_file(path.c_str(), "a"));
  REQUIRE_NOTHROW(assert_can_open_file(path.c_str(), "a+"));

  REQUIRE_THROWS_AS(assert_can_open_file(path.c_str(), "NONEXISTANT_FILE_MODE"),
                    std::runtime_error);

  // cleanup
  std::filesystem::remove_all(tmp);
}


TEST_CASE("Check char* equality") {
  REQUIRE(are_equal("A", "A"));
  REQUIRE(are_equal("ABCDEF", "ABCDEF"));
  REQUIRE(!are_equal("ABCDEF", "BACDEF"));
  REQUIRE(!are_equal("ABCDEF", "GHIJKL"));
}

TEST_CASE("tdms_vector_utils") {
  SECTION("has_elements") {
    vector<double> v_double;
    // Vector that is initialised with default behaviour should have 0 elements
    REQUIRE(!tdms_vector_utils::has_elements(v_double));
    // After resizing, the vector will now have elements
    v_double.resize(1);
    REQUIRE(tdms_vector_utils::has_elements(v_double));
  }

  SECTION("max") {
    vector<double> v_double = {1., 1.1, 1.2, 1.25, 1.2};
    REQUIRE(tdms_vector_utils::max(v_double) == Approx(1.25));
  }

  SECTION("index") {
    vector<int> v_int = {0, 1, 2, 1, 3};
    REQUIRE(tdms_vector_utils::index(v_int, 0) == 0);
    REQUIRE(tdms_vector_utils::index(v_int, 1) == 1);
    REQUIRE(tdms_vector_utils::index(v_int, 2) == 2);
    REQUIRE(tdms_vector_utils::index(v_int, 3) == 4);
    REQUIRE(tdms_vector_utils::index(v_int, 5) == -1);
  }

  SECTION("to_vector_int") {
    vector<double> v_double = {-1., 0., 1., 5.};
    vector<int> v_int = tdms_vector_utils::to_vector_int(v_double);
    bool converted_correctly =
            v_int[0] == -1 && v_int[1] == 0 && v_int[2] == 1 && v_int[3] == 5;
    REQUIRE(converted_correctly);
  }
}
