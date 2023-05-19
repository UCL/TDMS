#pragma once

#include <vector>

#include "cell_coordinate.h"
#include "matlabio.h"
#include "matrix.h"

/**
 * @brief Template class for storing two dimensional data as a strided vector.
 * @details Two dimensional data is stored as a strided vector, in row-major (C)
 * format. The column dimension has the fastest-varying index, so
 *
 * this->(i, j) = data_[i * n_cols_ + j].
 *
 * This is consistent with the way hdf5 files store array-like data, see
 * https://support.hdfgroup.org/HDF5/doc1.6/UG/12_Dataspaces.html.
 * @tparam T Numerical datatype
 */
template<typename T>
class Matrix {
  friend class HDF5Reader;
  friend class HDF5Writer;

protected:
  int n_rows_ = 0;
  int n_cols_ = 0;

  /*! Strided vector that will store the array data */
  std::vector<T> data_;

  /*! The total number of elements, according to the dimension values currently
   * set. */
  int n_elements() const { return n_rows_ * n_cols_; }

public:
  Matrix() = default;
  Matrix(int n_rows, int n_cols) { allocate(n_rows, n_cols); }

  int get_n_rows() const { return n_rows_; }
  int get_n_cols() const { return n_cols_; }

  /** @brief Subscript operator for the Matrix, retrieving the (i,j)-th element
   */
  T &operator()(int row, int col) { return data_[row * n_cols_ + col]; }
  /** @copydoc operator() */
  T operator()(int row, int col) const { return data_[row * n_cols_ + col]; }

  /**
   * @brief Allocate memory for this Matrix given the dimensions passed.
   *
   * @param n_rows,n_cols Number of rows and columns to assign.
   */
  void allocate(int n_rows, int n_cols) {
    n_rows_ = n_rows;
    n_cols_ = n_cols;
    data_.resize(n_elements());
  }

  /**
   * @brief Initialise this Matrix from a 2D-buffer of matching size.
   * @details Data values are copied so that membership over this->data_ is
   * preserved.
   * @param buffer 2D buffer to read from
   * @param n_rows,n_cols "Shape" of the read buffer to assign to this Matrix
   * @param buffer_leads_n_cols If true, the buffer to read from is assumed to
   * have dimensions [n_cols][n_rows]. If false, it is assumed to have
   * dimensions [n_rows][n_cols].
   */
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

  /** @brief Whether this Matrix contains any elements */
  bool has_elements() const { return n_elements() != 0; }
};

/** @brief Matrix of C-coefficients. See the pdf documentation for their
 * definition. */
typedef Matrix<double> CCoefficientMatrix;

/** TODO: Docstring */
struct GratingStructure : public Matrix<int> {
  GratingStructure(const mxArray *ptr, int I_tot);
};

/**
 * @brief Defines the numerical aperture of the objective, assuming that the
 * lens is centred on the origin of the PSTD simulation.
 * @details In particular, since the fibre modes are imaged onto a Fourier plane
 * of both the physical fibre and the sample, the field scattered by the sample
 * and collected by the objective lens can have only a finite spatial support in
 * the aperture of the objective lens.
 *
 * Pupil(i, j) thus takes the value 1 for those (i,j) indices within the
 * aperture of the lens.
 */
struct Pupil : public Matrix<double> {
  void initialise_from_matlab(const mxArray *ptr, int n_rows, int n_cols);
};

/**
 * @brief n_vertices-by-3 storage used for holding Yee cell indices.
 * @details Each row of a Vertices instance consists of three integers (i, j, k)
 * that form the index of a particular Yee cell in the simulation space.
 */
struct Vertices : public Matrix<int> {
  Vertices() = default;

  void initialise_from_matlab(const mxArray *ptr);

  int n_vertices() { return n_rows_; }

  /** @brief Convert the (i,j,k) index stored across the columns of row to an
   * ijk struct. */
  ijk index_in_row(int row) const;
};
