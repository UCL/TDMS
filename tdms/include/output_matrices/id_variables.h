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
 * The x and y members store the optical frequency (wavelength) dependent
 * detector signals that have been acquired in the Fourier plane of the
 * objective lens, corresponding to the respective components of the field.
 * These quantities are evaluated when exdetintegral is set to 1, and are
 * indexed by the wavenumbers and modes specified in the input file.
 *
 * The paper Peter R. T. Munro, "Exploiting data redundancy in computational
 * optical imaging," Opt. Express 23, 30603-30617 (2015), explains the
 * underlying theory and. In particular, how the coupling coefficient of a field
 * coupled into an optical fiber can be evaluated at different planes within the
 * optical system.
 *
 * Whilst a single mode optical fiber has just one mode which the field must
 * match in order to be coupled into the fibre, this functionality allows the
 * user to specify a multi-mode fibre, or a single-mode fiber that can be
 * translated to different positions.
 *
 * The ID output is a MATLAB struct that we create during runtime, so we need to
 * free the memory we reserve upon destruction, as well as link pointers to the
 * MATLAB structs to our native C++ arrays.
 */
class IDVariables {
private:
  int n_frequencies = 0;//< Number of frequencies that we're extracting at
  int n_det_modes = 0;  //< The number of fibre modes

  bool memory_assigned = false;//< Flags whether MATLAB memory has been assigned
                               // and needs to be free'd
public:
  mxArray *x_ptr =
          nullptr;//< Holds the array in the Idx field of OutputMatrices["Id"]
  mxArray *y_ptr =
          nullptr;//< Holds the arrays in the Idy field of OutputMatrices["Id"]

  /*! Indexed by (i_wavenum, i_m).
  i_wavenum indexes the wavenumbers specified in the input file using the
  k_vec.i_m indices in the main loop. i_m indexes the different modes specified
  by the input detmodevec. */
  std::complex<double> **x = nullptr, **y = nullptr;

  // pointers to the real (re) and imaginary (im) parts of the data in the Id
  // output
  double **x_real = nullptr, **y_real = nullptr, **x_imag = nullptr,
         **y_imag = nullptr;

  IDVariables() = default;

  /**
   * @brief Creates MATLAB arrays that constitute the fields of the ID output,
   * and assigns them as the appropriate fields of the ID output.
   *
   * The ID output is a 1-by-1 structure array with 2 fields, which is pointed
   * to in id_pointer. In order to assign arrays as the field values, we need to
   * create fresh MATLAB structures using the member variables of this class.
   * The arrays assigned to the two fieldnames, "x" and "y", are of side
   * n_frequencies-by-n_det_modes.
   *
   * @param id_pointer Pointer to the ID output structure array
   * @param _n_frequencies The number of frequencies we are considering in this
   * simulation
   * @param n_det_modes Number of fibre modes specified
   */
  void link_to_pointer(mxArray *&id_pointer, int _n_frequencies,
                       int n_det_modes);

  ~IDVariables();
};
