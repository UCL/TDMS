/**
 * @file eh_vec.h
 * @author William Graham (william.graham@ucl.ac.uk)
 * @brief Declaration for the storage class used when performing Fourier
 * transforms via fftw
 */
#pragma once

#include <fftw3.h>

/**
 * @brief Storage class for Fourier transform values, obtained when performing
 * the pseudo-spectral timestep.
 * @details Storage is organised as a 2D buffer of fftw_complex types, of size
 * n_rows_ by n_cols. Elements can be accessed via the [] operator.
 */
class EHVec {
protected:
  int n_rows_ = 0;
  int n_cols_ = 0;

  fftw_complex **data_ = nullptr;

  /** @brief Free the memory assigned to this instance */
  void free_memory();

public:
  EHVec() = default;

  fftw_complex *operator[](int col_index) { return data_[col_index]; }

  /**
   * @brief Allocate storage for this instance. Existing storage will be lost
   * (if present).
   * @details Reallocation of an already-allocated EHVec frees the current
   * buffer before reassigning the new buffer. Any values stored in the old
   * buffer are lost, EVEN IF THE NEW BUFFER IS OF GREATER SIZE THAN THE OLD
   * BUFFER.
   * @param n_rows,n_cols Shape of storage to allocate.
   */
  void allocate(int n_rows, int n_cols);

  /** @brief Determine whether the instance has memory allocated. */
  bool is_allocated() const { return data_ != nullptr; }

  ~EHVec() { free_memory(); }
};
