/**
 * @file fdtd_grid_initialiser.h
 * @brief Initialisation of the FDTD grid
 */
#pragma once

#include <string>
#include <vector>

#include "mat_io.h"

/**
 * @brief A class to initialise an FDTD grid from a MATLAB file
 */
class fdtdGridInitialiser {

private:
  const mxArray *pointer = nullptr;          //< Pointer to the array
  const char *mat_filename = nullptr;        //< Filename of the MATLAB file
  std::vector<mwSize> dimensions = {0, 0, 0};//< The dimensions of the array

  /**
   * @brief Get a value from a integer attribute of the FDTD grid defined in a
   * .mat file
   *
   * @param key the name of the attribute
   * @return int the value of the attribute
   */
  mwSize value_of_attribute(const std::string &key);

public:
  /** @brief Construct a new fdtd Grid Initialiser object */
  fdtdGridInitialiser() {}
  /**
   * @brief Construct a new fdtd Grid Initialiser object
   *
   * @param fdtd_pointer pointer to the FDTD grid
   * @param mat_filename the filename of the MATLAB file
   */
  fdtdGridInitialiser(const mxArray *fdtd_pointer, const char *mat_filename);

  /**
   * @brief Set an FDTD grid attribute to a tensor full of zeros
   * @param name of the attribute
   */
  void add_tensor(const std::string &name);
};
