/**
 * @file id_variables.h
 * @brief Contains the IDVariables class, handling the ID output matrix.
 */
#pragma once

#include <complex>

#include "matrix.h"

/**
 * @brief Class that handles the variables associated with the ID output.
 *
 * The ID output is a MATLAB struct that we create during runtime, so we need to free the memory we reserve upon destruction, as well as link pointers to the MATLAB structs to our native C++ arrays.
 */
class IDVariables {
private:
  int n_frequencies = 0;//< Number of frequencies that we're extracting at. ( = to f_ex_vec.size() )
  int n_det_modes = 0;//< D_tilde.num_det_modes()

  bool memory_assigned = false;//< Flags whether MATLAB memory has been assigned and needs to be free'd
public:
  mxArray *x_ptr = nullptr;//< Holds the array in the Idx field of OutputMatrices["Id"]
  mxArray *y_ptr = nullptr;//< Holds the arrays in the Idy field of OutputMatrices["Id"]

  std::complex<double> **x = nullptr, **y = nullptr;

  // pointers to the real (re) and imaginary (im) parts of the data in the Id output
  double **x_real = nullptr, **y_real = nullptr, **x_imag = nullptr, **y_imag = nullptr;

  IDVariables() = default;

  void link_to_pointer(mxArray *&id_pointer, int _n_frequencies, int n_det_modes);

  ~IDVariables();
};
