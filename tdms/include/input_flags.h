/**
 * @file input_flags.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Organises enumerated constants, names, and classes for handling
 * flag-variables that are passed to TDMS in the input file. Such variables do
 * not contain numerical data for use in the simulation, but affect the manner
 * in which the simulation is performed.
 */
#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>

#include "mat_io.h"

//! The maximum number of flags that could be present in the input file
#define NFLAGS 2

//! Namespace encompassing variables and enums related to flag values to be read
//! from the input file. Enums classes are set so that the value corresponding
//! to "true" is the corresponding option that is used when that flag is passed
//! in the input file.
namespace tdms_flags {
//! Lists the flag-variables that can be present in the input file, but are not
//! required to be present
const std::vector<std::string> flag_variables = {"use_pstd", "use_bli"};

//! The timestepping method to be used to forward-propagate the simulation.
//! use_pstd is the corresponding flag.
enum SolverMethod : bool { PseudoSpectral = true, FiniteDifference = false };
//! The interpolation method to use when extracting field values at Yee
//! cell centres. use_bli is the corresponding flag.
enum InterpolationMethod : bool { BandLimited = true, Cubic = false };
}// namespace tdms_flags

class InputFlags {
private:
  //! Name of the input file to read flags from
  std::string input_filename = "input_file.mat";
  //! Flag values present in the input file
  bool flag_values[NFLAGS] = {false};

  /**
   * @brief Returns the internal index in flag_values of the flag with the name
   * provided.
   *
   * @param flag_name Name of the flag
   * @return const int Index of this flag's value in flag_values
   */
  const int position_from_name(const std::string &flag_name) const {
    auto position = std::find(tdms_flags::flag_variables.begin(),
                              tdms_flags::flag_variables.end(), flag_name);
    if (position == tdms_flags::flag_variables.end()) {
      // Could not find the flag in the list of expected flags, return an error
      throw std::runtime_error(flag_name + " is not an expected flag");
    }
    return std::distance(tdms_flags::flag_variables.begin(), position);
  }

  /**
   * @brief Fetch the value of the flag variable from the input file.
   *
   * @param flag_name The name of the flag to fetch
   * @param fail_on_not_found If true, throw an error if the flag cannot be
   * located. Otherwise, return the default flag value.
   * @return true This flag was present and set to true
   * @return false This flag was present and set to false, OR was not present
   */
  bool fetch_flag_value(const std::string flag_name,
                        bool fail_on_not_found = false) const;

public:
  InputFlags() = default;
  /**
   * @brief Construct a new InputFlags object by reading flags from the input
   * file provided.
   *
   * The user may optionally specify that an error should be thrown if flags
   * cannot be found through the fail_on_not_found argument.
   *
   * @param _input_filename The input file to read from
   * @param fail_on_not_found If true, any flags that cannot be located in the
   * input file will cause an error to be thrown.
   */
  InputFlags(const std::string &_input_filename, bool fail_on_not_found = false)
      : input_filename(_input_filename) {
    for (const std::string &flag_name : tdms_flags::flag_variables) {
      flag_values[position_from_name(flag_name)] =
              fetch_flag_value(flag_name, fail_on_not_found);
    }
  };

  /**
   * @brief Return the value of the flag corresponding to the name provided.
   */
  bool operator[](const std::string &flag_name) const {
    return flag_values[position_from_name(flag_name)];
  }

  /** @brief Prompts spdlog to report the value of every expected flag in the
   * input file */
  void report_flag_state() {
    for (const std::string &flag : tdms_flags::flag_variables) {
      spdlog::info("Read " + flag + ": {}",
                   flag_values[position_from_name(flag)]);
    }
  }
};
