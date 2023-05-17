/**
 * @file test_argument_parser.cpp
 * @brief Tests of the argument parsing and file I/O.
 */
#include "argument_parser.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include <catch2/catch_test_macros.hpp>

using std::string;
using std::vector;

/**
 * @brief Test that the argument namespace recovers the input arguments
 * provided.
 */
TEST_CASE("Test namespace") {

  SECTION("Help") {
    const char *input_args[] = {"tdms", "-h"};
    auto args = ArgumentNamespace(2, const_cast<char **>(input_args));
    REQUIRE(args.have_flag("-h"));
  }
  SECTION("Version") {
    const char *input_args[] = {"tdms", "--version"};
    auto args = ArgumentNamespace(2, const_cast<char **>(input_args));
    REQUIRE(args.have_flag("--version"));
  }
}

/** Drop the first and last component of an std::vector (needs size > 2!!). */
template<typename T>
vector<T> drop_first_and_last(const vector<T> &v) {
  return vector<T>(v.begin() + 1, v.end() - 1);
}

/** Convert a vector of std::strings to an array of chars */
const char **vector_to_array(const vector<string> &vec) {
  const char **arr = new const char *[vec.size()] {};
  std::transform(vec.begin(), vec.end(), arr,
                 [](const std::string &s) { return s.c_str(); });
  return arr;
}
// TODO: check with Will that this doesn't shadow the implementation

/**
 * @brief Test that the argument parser returns the correct filenames.
 */
TEST_CASE("Test parsing two filenames") {

  const vector<string> input_args = {"tdms", "input_file.mat",
                                     "output_file.mat"};
  ArgumentParser ap;
  auto args = ap.parse_args(input_args.size(),
                            const_cast<char **>(vector_to_array(input_args)));

  SECTION("Inputs") {
    vector<string> expected_filenames = drop_first_and_last(input_args);
    vector<string> parsed_filenames = args.input_filenames();
    REQUIRE(expected_filenames == parsed_filenames);
  }
  SECTION("Output") {
    string parsed = args.output_filename();
    string expected = input_args[input_args.size() - 1];
    REQUIRE(expected == parsed);
  }
  SECTION("Grid") {
    REQUIRE(!args.has_grid_filename());
    REQUIRE_THROWS_AS(args.grid_filename(), std::runtime_error);
  }
}

/**
 * @brief Test that the argument parser returns the correct filenames.
 */
TEST_CASE("Test parsing three filenames") {
  vector<string> input_args = {"tdms", "input_file.mat", "grid_file.mat",
                               "output_file.mat"};
  ArgumentParser ap;
  auto args = ap.parse_args(input_args.size(),
                            const_cast<char **>(vector_to_array(input_args)));

  SECTION("Input") {
    vector<string> expected_filenames = drop_first_and_last(input_args);
    vector<string> parsed_filenames = args.input_filenames();
    REQUIRE(expected_filenames == parsed_filenames);
  }
  SECTION("Output") {
    string parsed = args.output_filename();
    string expected = input_args[input_args.size() - 1];
    REQUIRE(expected == parsed);
  }
  SECTION("Grid") {
    string parsed = args.grid_filename();
    string expected = input_args[input_args.size() - 2];
    REQUIRE(args.has_grid_filename());
    REQUIRE(expected == parsed);
  }
}

/**
 * @brief Test that the argument correctly errors with the wrong number of file
 * names.
 */
TEST_CASE("Test wrong number of file names") {
  SECTION("Zero") {
    const char *input_args[] = {"tdms"};
    auto args = ArgumentNamespace(1, const_cast<char **>(input_args));
    REQUIRE(!args.have_correct_number_of_filenames());
    REQUIRE_THROWS_AS(args.input_filename(), std::runtime_error);
    REQUIRE_THROWS_AS(args.grid_filename(), std::runtime_error);
    REQUIRE_THROWS_AS(args.output_filename(), std::runtime_error);
  }
  SECTION("One") {
    const char *input_args[] = {"tdms", "input_file.mat"};
    auto args = ArgumentNamespace(2, const_cast<char **>(input_args));
    REQUIRE(!args.have_correct_number_of_filenames());
    REQUIRE_THROWS_AS(args.grid_filename(), std::runtime_error);
    REQUIRE_THROWS_AS(args.output_filename(), std::runtime_error);
  }
  SECTION("Four") {
    const char *input_args[] = {"tdms", "input_file.mat", "grid_file.mat",
                                "output_file.mat", "fourth_file.mat"};
    auto args = ArgumentNamespace(5, const_cast<char **>(input_args));
    REQUIRE(!args.have_correct_number_of_filenames());
  }
}

/** @brief Test of the CLI options. */
TEST_CASE("CL options") {
  vector<string> input_args = {"tdms", "input_file.mat", "output_file.mat"};
  SECTION("Doesn't have flags") {
    auto args =
            ArgumentNamespace(input_args.size(),
                              const_cast<char **>(vector_to_array(input_args)));
    // ArgumentNamespace should ignore executable name when counting input args
    REQUIRE(args.num_non_flag == 2);
  }
}
