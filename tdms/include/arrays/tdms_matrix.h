#pragma once

#include <vector>

#include "matlabio.h"
#include "matrix.h"

template<typename T>
class TDMSMatrix {
protected:
  int n_rows_ = 0;
  int n_cols_ = 0;

  std::vector<T> data_;

  int n_elements() const { return n_rows_ * n_cols_; }

public:
  TDMSMatrix() = default;
  TDMSMatrix(int n_rows, int n_cols) { allocate(n_rows, n_cols); }

  T &operator()(int row, int col) {
    return data_[row * n_cols_ + col];
  }// TODO: convert old syntax :(

  void allocate(int n_rows, int n_cols) {
    n_rows_ = n_rows;
    n_cols_ = n_cols;

    data_.resize(n_elements());
  }

  void initialise(T **buffer, int n_rows, int n_cols,
                  bool buffer_leads_n_cols = false) {
    allocate(n_rows, n_cols);
    if (buffer_leads_n_cols) {
      for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) { operator()(i, j) = buffer[j][i]; }
      }
    } else {
      for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) { operator()(i, j) = buffer[i][j]; }
      }
    }
  }

  bool has_elements() const { return n_elements() != 0; }
};

class CCoefficientMatrix : public TDMSMatrix<double> {};

struct GratingStructure : public TDMSMatrix<int> {
  GratingStructure(const mxArray *ptr, int I_tot);
};

// class Pupil : public Matrix<double> {
// public:
//   Pupil() = default;

//   void initialise(const mxArray *ptr, int n_rows, int n_cols);

//   ~Pupil();
// };

// class Vertices : public Matrix<int> {
// public:
//   Vertices() = default;

//   void initialise(const mxArray *ptr);

//   int n_vertices() { return n_rows_; }

//   ~Vertices() {
//     if (has_elements()) { free_cast_matlab_2D_array(matrix); }
//     matrix = nullptr;
//   };
// };
