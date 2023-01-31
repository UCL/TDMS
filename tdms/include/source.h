/**
 * @file source.h
 */
#pragma once

#include <string>

#include "mat_io.h"

class Source {
private:
  bool no_data_stored = true;//!< Flags if the array is empty to avoid pointer preservation
public:
  double ***real = nullptr;//!< Real data for the source term
  double ***imag = nullptr;//!< Imag data for the source term

  Source(const mxArray *ptr, int dim1, int dim2, const std::string &name);

  /** @brief Check if the source term is empty (true) or not (false) */
  bool is_empty() { return no_data_stored; }
};
