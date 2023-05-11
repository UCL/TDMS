/**
 * @file dtilde.h
 * @brief Container class for the DTilde object, storing the fibre modes in the
 * Fourier plane of the objective lens.
 */
#pragma once

#include <complex.h>
#include <string>
#include <vector>

#include "arrays/tensor3d.h"
#include "matrix.h"

/**
 * @brief Stores the fibre modes in the Fourier plane of the objective lens.
 * @details The "Tilde" indicates that these quantities are in a Fourier plane
 * relative to where the optical fibre is actually located, meaning that it has
 * a Fourier relationship relative to the physical fibre mode(s).
 */
class DTilde {
protected:
  int n_det_modes = 0;//<! Number of modes specified

  /**
   * @brief Set one of the two components by reading from an existing MATLAB
   * buffer.
   *
   * @param tensor Component to read into
   * @param ptr Pointer to the [struct] memory buffer to read from
   * @param name [Debug] Name of the field to read from
   * @param n_rows,n_cols Dimensions to assign to the component
   */
  static void set_component(Tensor3D<std::complex<double>> &tensor,
                            const mxArray *ptr, const std::string &name,
                            int n_rows, int n_cols);

public:
  DTilde() = default;

  /** @brief (Fetch the) number of fibre modes specified */
  int num_det_modes() const { return n_det_modes; };

  Tensor3D<std::complex<double>> x;
  Tensor3D<std::complex<double>> y;

  /**
   * @brief Initialise the fibre modes by reading from a MATLAB buffer.
   *
   * @param ptr Pointer to a MATLAB struct containing the data to read
   * @param n_rows,n_cols Dimensions to assign to the component
   */
  void initialise(const mxArray *ptr, int n_rows, int n_cols);
};
