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
 * @details Three dimensional data is stored as a strided vector, in row-major
 * (C) format. The last listed dimension has the fastest varying index, and the
 * first listed dimension has the slowest listed. Explicitly,
 *
 * this->(i, j, k) = this[i*n_cols*n_layers + j*n_layers + k].
 *
 * This is consistent with the way hdf5 files store array-like data, see
 * https://support.hdfgroup.org/HDF5/doc1.6/UG/12_Dataspaces.html.
 * @tparam T Numerical datatype
 */
template<typename T>
class Tensor3D {
private:
  /** @brief Convert a 3D (i,j,k) index to the corresponding index in the
   * strided storage. */
  int to_global_index(int i, int j, int k) const {
    return i * n_layers_ * n_cols_ + j * n_layers_ + k;
  }
  int to_global_index(const ijk &index_3d) const {
    return to_global_index(index_3d.i, index_3d.j, index_3d.k);
  }

protected:
  int n_layers_ = 0;
  int n_cols_ = 0;
  int n_rows_ = 0;

  /*! Strided vector that will store the array data */
  std::vector<T> data_;

  /*! The total number of elements, according to the dimension values currently
   * set. */
  int total_elements() const { return n_layers_ * n_cols_ * n_rows_; }

public:
  Tensor3D() = default;
  Tensor3D(int n_layers, int n_cols, int n_rows) {
    allocate(n_layers, n_cols, n_rows);
  }
  Tensor3D(T ***buffer, int n_layers, int n_cols, int n_rows) {
    initialise(buffer, n_layers, n_cols, n_rows);
  }

  /** 
   * @brief Subscript operator for the Tensor, retrieves the (i,j,k)-th element.
   * @note Does not 'call' the tensor, but rather accesses the data.
   * @example ... 
   * @seealso Tensor3D::operator[](const ijk&)
   */
  T &operator()(int i, int j, int k) { return data_[to_global_index(i, j, k)]; }
  T operator()(int i, int j, int k) const {
    return data_[to_global_index(i, j, k)];
  }

  /** @brief Subscript operator for the Tensor, retrieving the (i,j,k)-th
   * element. */
  T &operator[](const ijk &index_3d) {
    return data_[to_global_index(index_3d)];
  }
  T operator[](const ijk &index_3d) const {
    return data_[to_global_index(index_3d)];
  }

  /**
   * @brief Whether the tensor contains any elements.
   * @details Strictly speaking, checks whether the total_elements() method
   * returns 0, indicating that this tensor has no size and thus no elements.
   */
  bool has_elements() const { return total_elements() != 0; }

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
    data_.resize(total_elements());
  }

  /**
   * @brief Initialise this tensor from a 3D-buffer of matching size.
   * @details Data values are copied so that membership over this->data() is
   * preserved.
   * @param buffer 3D buffer to read from
   * @param n_layers,n_cols,n_rows "Shape" of the read buffer to assign to this
   * tensor.
   * @param buffer_leads_n_layers If true, the buffer to read from is
   * assumed to have dimensions [n_layers][n_cols][n_rows]. If false, it is
   * assumed to have the dimensions in reverse order,
   * [n_rows][n_cols][n_layers].
   */
  void initialise(T ***buffer, int n_layers, int n_cols, int n_rows,
                  bool buffer_leads_n_layers = false) {
    allocate(n_layers, n_cols, n_rows);
    if (buffer_leads_n_layers) {
      for (int k = 0; k < n_layers_; k++) {
        for (int j = 0; j < n_cols_; j++) {
          for (int i = 0; i < n_rows_; i++) {
            operator()(i, j, k) = buffer[k][j][i];
          }
        }
      }
    } else {
      for (int k = 0; k < n_layers_; k++) {
        for (int j = 0; j < n_cols_; j++) {
          for (int i = 0; i < n_rows_; i++) {
            operator()(i, j, k) = buffer[i][j][k];
          }
        }
      }
    }
  }

  /** @brief Set all elements in the tensor to 0. */
  void zero() { std::fill(data_.begin(), data_.end(), 0); }

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
    for (const T &element_value : data_) {
      norm_val += std::abs(element_value) * std::abs(element_value);
    }
    return std::sqrt(norm_val);
  }
};
