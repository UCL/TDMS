/**
 * @file argument_parser.h
 * @brief Parse the command line options.
 */
#pragma once

#include <string>
#include <vector>

#include "globals.h"

/**
 * @brief Wraps a vector of string CL arguments with some helpful functionality.
 */
class ArgumentNamespace {

private:
  std::vector<std::string> arguments;         //< The arguments from flags
  std::vector<std::string> non_flag_arguments;//< Arguments not from flags

public:
  int num_non_flag = 0;//< A count of arguments not from flags

  /**
   * @brief Construct a new Argument Namespace object
   *
   * @param n_args the number of arguments (argc)
   * @param argv pointers to the arguments (argv)
   */
  explicit ArgumentNamespace(int n_args, char *argv[]);

  /**
   * @brief Searches the arguments for the flag provided
   *
   * @param flag to search for
   * @return true if the flag is present
   * @return false otherwise
   */
  bool have_flag(std::string const &flag) const;

  /**
   * @brief Return true if the user has requested the output file we written in
   * compressed format.
   * @details Compressed output is toggled by the -m, --minimise-file-size
   * options on the CLI.
   *
   * This option is used when the outputs are being written, after the
   * simulation has been performed.
   */
  bool compressed_output() const;

  /**
   * @brief Have we been provided with a grid filename?
   *
   * @return true if there are 3 non-flag arguments (therefore a grid filename)
   * @return false otherwise
   */
  bool has_grid_filename() const;

  /**
   * @brief Check that the correct number of filename arguments are provided
   * (either 2 or 3 non-flag arguments)
   *
   * @return true if correct
   * @return false otherwise
   */
  bool have_correct_number_of_filenames() const;

  /**
   * @brief Check whether an argument is a flag (starts '-')
   *
   * @param arg the argument
   * @return true if arg is a flag
   * @return false otherwise
   */
  static bool is_a_flag_argument(std::string arg);

  /**
   * @brief Gets the input filename
   *
   * @return const char* the input filename
   */
  const char *input_filename();

  /**
   * @brief Gets the output filename
   *
   * The output filename is either the second or third positional argument
   * depending on whether a grid filename is provided or not.
   *
   * @return const char*  the output filename
   */
  const char *output_filename();

  /**
   * @brief Gets the grid filename
   *
   * @return const char* the grid filename
   */
  const char *grid_filename();

  /**
   * @brief Get all input filenames
   *
   * A vector containing the input filename and the grid filename (if
   * provided)
   *
   * @return std::vector<std::string> the input filenames
   */
  std::vector<std::string> input_filenames();

  /**
   * @brief Check that all input and output files can be accessed with the
   * correct privileges
   */
  void check_files_can_be_accessed();
};

/**
 * @brief Performs the argument parsing and returns an ArgumentNamespace
 */
class ArgumentParser {

private:
  /** Prints the help message (all options).  */
  static void print_help_message();
  /** Prints the version of the compiled executable. */
  static void print_version();

public:
  /**
   * @brief Parse the command line arguments and perform relevant actions.
   *
   * @param n_args The number of arguments (argc)
   * @param arg_ptrs Pointers to to the arguments (argv)
   * @return ArgumentNamespace populated with options
   */
  static ArgumentNamespace parse_args(int n_args, char *arg_ptrs[]);
};
