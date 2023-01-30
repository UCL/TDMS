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

  /**
   * @brief Creates MATLAB arrays that constitute the fields of the ID output, and assigns them as the appropriate fields of the ID output.
   *
   * The ID output is a 1-by-1 structure array with 2 fields, which is pointed to in id_pointer.
   * In order to assign arrays as the field values, we need to create fresh MATLAB structures using the member variables of this class.
   * The arrays assigned to the two fieldnames, "x" and "y", are of side n_frequencies-by-n_det_modes.
   *
   * @param id_pointer Pointer to the ID output structure array
   * @param _n_frequencies The number of frequencies we are considering in this simulation
   * @param n_det_modes TODO:: DTilde.n_det_modes()
   */
  void link_to_pointer(mxArray *&id_pointer, int _n_frequencies, int n_det_modes);

  ~IDVariables();
};
