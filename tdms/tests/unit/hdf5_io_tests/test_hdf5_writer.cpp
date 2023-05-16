/**
 * @file test_hdf5_writer.cpp
 * @author Sam Cunliffe, William Graham
 * @brief Unit tests for HDF5Writer output functions.
 */
#include "hdf5_io/hdf5_reader.h"
#include "hdf5_io/hdf5_writer.h"

#include <filesystem>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <spdlog/spdlog.h>

#include "unit_test_utils.h"

using Catch::Approx;
using namespace std;
using tdms_tests::create_tmp_dir;

TEST_CASE("HDF5Writer: Write doubles to a group.") {

  // setup
  string tmp_dir = create_tmp_dir();
  string file_name = tmp_dir + "/some_test_data.hdf5";
  spdlog::info("Temporary file: {}", file_name);

  string group_name = "group_name";
  string double_dataset = "array_of_10_doubles";
  string int_dataset = "array_of_5_ints";

  {
    HDF5Writer writer(file_name);
    {
      vector<double> write_me_out = {.0, .1, .2, .3, .4, .5, .6, .7, .8, .9};
      writer.write_dataset_to_group(group_name, double_dataset, write_me_out);
    }
    {
      vector<int> write_me_out = {0, 1, 2, 3, 4, 5};
      writer.write_dataset_to_group(group_name, int_dataset, write_me_out);
    }
    REQUIRE(writer.contains(group_name));
  }

  // Read data back to confirm entries are as expected
  {
    HDF5Reader reader(file_name);
    REQUIRE(reader.contains(group_name));
    {
      vector<double> read_back_in;
      reader.read_dataset_in_group(group_name, double_dataset, read_back_in);

      bool all_entries_correct = true;
      for (int i = 0; i < 10; i++) {
        all_entries_correct = all_entries_correct &&
                              (read_back_in[i] == Approx((double) i / 10.));
      }
      REQUIRE(all_entries_correct);
    }
    // TODO: fix this!
    // {
    //  vector<int> read_back_in;
    //  reader.read_dataset_in_group(group_name, double_dataset, read_back_in);

    //  bool all_entries_correct = true;
    //  for (int i = 0; i < 5; i++) {
    //    all_entries_correct = all_entries_correct &&
    //                          (read_back_in[i] == i);
    //  }
    //  REQUIRE(all_entries_correct);
    // }
  }

  // Teardown
  filesystem::remove_all(tmp_dir);
}
