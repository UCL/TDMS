/**
 * @file tensor3d.h
 * @author William Graham
 * @brief template class for storing three-dimensional data arrays.
 */
#pragma once

#include <vector>

#include "cell_coordinate.h"

/**
 * @brief Template class for storing three-dimensional data as a strided vector.
 * @details Three dimensional data is stored as a strided vector, with the
 * layers index being the slowest varying. That is, given the number of layers,
 * columns, and rows as n_layers, n_cols, n_rows respectively:
 *
 * this->(i, j, k) = this[k*n_columns*n_rows + j*n_rows + i].
 *
 * This is consistent with the format in which array-like data is stored in hdf5
 * files, and MATLAB files.
 *
 * @tparam T Numerical datatype
 */
template<typename T>
class Tensor3D : public std::vector<T> {
protected:
  int n_layers_ = 0;
  int n_cols_ = 0;
  int n_rows_ = 0;

  /*! The total number of elements, according to the dimension values currently
   * set. */
  int total_elements() const { return n_layers_ * n_cols_ * n_rows_; }

public:
  Tensor3D() : std::vector<T>(){};
  Tensor3D(int n_layers, int n_cols, int n_rows) {
    allocate(n_layers, n_cols, n_rows);
  }
  Tensor3D(T ***buffer, int n_layers, int n_cols, int n_rows) {
    initialise(buffer, n_layers, n_cols, n_rows);
  }

  /** @brief Subscript operator for the Tensor, retrieving the (i,j,k)-th
   * element. */
  T &operator()(int i, int j, int k) {
    return *(this->begin() + k * n_cols_ * n_rows_ + j * n_rows_ + i);
  }
  T operator()(int i, int j, int k) const {
    return *(this->begin() + k * n_cols_ * n_rows_ + j * n_rows_ + i);
  }

  /** @brief Subscript operator for the Tensor, retrieving the (i,j,k)-th
   * element. */
  T &operator[](const ijk &index_3d) {
    return this->operator()(index_3d.i, index_3d.j, index_3d.k);
  }
  T operator[](const ijk &index_3d) const {
    return this->operator()(index_3d.i, index_3d.j, index_3d.k);
  }

  /**
   * @brief Whether the tensor contains any elements.
   * @details Strictly speaking, checks whether the total_elements() method
   * returns 0, indicating that this tensor has no size and thus no elements.
   */
  bool has_elements() const { return this->total_elements() != 0; }

  /**
   * @brief Allocate memory for this tensor given the dimensions passed.
   *
   * @param n_layers,n_cols,n_rows Number of layers, columns, and rows to
   * assign.
   */
  void allocate(int n_layers, int n_cols, int n_rows) {
    n_layers_ = n_layers;
    n_cols_ = n_cols;
    n_rows_ = n_rows;
    this->reserve(this->total_elements());
  }

  /**
   * @brief Initialise this tensor from a 3D-buffer of matching size.
   * @details Data values are copied so that membership over this->data() is
   * preserved.
   * @param buffer 3D buffer to read from
   * @param n_layers,n_cols,n_rows "Shape" of the read buffer to assign to this
   * tensor.
   */
  void initialise(T ***buffer, int n_layers, int n_cols, int n_rows) {
    this->allocate(n_layers, n_cols, n_rows);
    for (int k = 0; k < n_layers_; k++) {
      for (int j = 0; j < n_cols_; j++) {
        for (int i = 0; i < n_rows_; i++) {
          this->operator()(i, j, k) = buffer[k][j][i];
        }
      }
    }
  }

  /** @brief Set all elements in the tensor to 0. */
  void zero() { std::fill(this->begin(), this->end(), 0); }

  /**
   * @brief Computes the Frobenius norm of the tensor.
   * @details The Frobenius norm is defined as:
   *
   * fro_norm = \f$\sqrt{
   * \sum_{i=0}^{n_rows}\sum_{j=0}^{n_cols}\sum_{k=0}^{n_layers} |t[k][j][i]|^2
   * }\f$.
   *
   * In practice, our data is flat and since addition is commutative, we can
   * just iterate over the flattened data structure instead of using nested for
   * loops.
   */
  double frobenius() const {
    T norm_val = 0;
    for (const T &element_value : *this) {
      norm_val += abs(element_value) * abs(element_value);
    }
    return sqrt(norm_val);
  }
};
