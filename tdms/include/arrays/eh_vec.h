#pragma once

#include <fftw3.h>
#include <spdlog/spdlog.h>

// EHVec[n_cols][n_rows]
class EHVec {
protected:
  int n_rows_ = 0;
  int n_cols_ = 0;

  fftw_complex **data_ = nullptr;

  void free_memory() {
    spdlog::info("Freeing memory");
    if (is_allocated()) {
      for (int i = 0; i < n_rows_; i++) { std::free(data_[i]); }
    }
    std::free(data_);
    data_ = nullptr;
    spdlog::info("Free'd");
  }

public:
  EHVec() = default;

  fftw_complex *operator[](int col_index) { return data_[col_index]; }

  void allocate(int n_rows, int n_cols) {
    // If already allocated we need to clean up old memory before resizing!
    if (is_allocated()) {
      spdlog::warn("EHVec is already assigned - freeing old memory before "
                   "reallocating!");
      free_memory();
    }

    n_rows_ = n_rows;
    n_cols_ = n_cols;

    data_ = (fftw_complex **) malloc(n_rows_ * sizeof(fftw_complex *));
    for (int i = 0; i < n_cols_; i++) {
      data_[i] = (fftw_complex *) malloc(n_cols_ * sizeof(fftw_complex));
    }
  }

  bool is_allocated() const { return data_ != nullptr; }

  ~EHVec() { free_memory(); }
};
