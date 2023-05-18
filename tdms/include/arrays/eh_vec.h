#pragma once

#include <fftw3.h>
#include <spdlog/spdlog.h>

// EHVec[n_cols][n_rows]
class EHVec {
protected:
  int n_rows_ = 0;
  int n_cols_ = 0;

  fftw_complex **data_ = nullptr;

  void free_memory();

public:
  EHVec() = default;

  fftw_complex *operator[](int col_index) { return data_[col_index]; }

  void allocate(int n_rows, int n_cols);

  bool is_allocated() const { return data_ != nullptr; }

  ~EHVec() { free_memory(); }
};
