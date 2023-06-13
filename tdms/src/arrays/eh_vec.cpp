#include "arrays/eh_vec.h"

#include <spdlog/spdlog.h>

void EHVec::free_memory() {
  if (is_allocated()) {
    for (int i = 0; i < n_rows_; i++) { fftw_free(data_[i]); }
    free(data_);
  }
  data_ = nullptr;
}

void EHVec::allocate(int n_rows, int n_cols) {
  // If already allocated we need to clean up old memory before resizing!
  if (is_allocated()) {
    spdlog::warn("EHVec is already assigned - freeing old memory before "
                 "reallocating!");
    free_memory();
  }

  n_rows_ = n_rows;
  n_cols_ = n_cols;

  data_ = (fftw_complex **) malloc(n_rows_ * sizeof(fftw_complex *));
  for (int i = 0; i < n_rows_; i++) {
    data_[i] = (fftw_complex *) malloc(n_cols_ * sizeof(fftw_complex));
  }
}
