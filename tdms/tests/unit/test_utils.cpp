/**
 * @file test_utils.cpp
 * @brief Tests of the utility functions
 */
#include "utils.h"

#include <fstream>
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include "unit_test_utils.h"

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

TEST_CASE("Test threeway max") {
  REQUIRE(max(1, 2, 3) == 3);
  REQUIRE(max(3, 2, 1) == 3);
  REQUIRE(max(2, 3, 1) == 3);
  REQUIRE(max(2, 1, 3) == 3);
  REQUIRE(max(2, 2, 3) == 3);
  REQUIRE(max(3, 3, 3) == 3);
  REQUIRE(max(1, 100, 10) == 100);
}
