/**
 * @file dimensions.h
 * @brief
 */
#pragma once

#include "mat_io.h"

/**
 * @brief A class encapsulating the dimensions in i, j, k.
 */
class Dimensions {
private:
  int i = 0;//!< The size (number of elements) in the x direction.
  int j = 0;//!< The size (number of elements) in the y direction.
  int k = 0;//!< The size (number of elements) in the z direction.

  /**
   * @brief Is this describing an nD object?
   *
   * @param n The number of dimensions to test.
   * @return true If the object is n-dimensional.
   * @return false Otherwise.
   */
  bool are_nd(int n) const { return (bool(i) + bool(j) + bool(k)) == n; }

public:
  /**
   * @brief Access the size in each dimension 0=i, 1=j, 2=k.
   *
   * @param value The index.
   * @return int The size.
   */
  int operator[](int value) const;

  /**
   * @brief Construct a new Dimensions object from a pointer to a MATLAB array.
   *
   * @param ptr The MATLAB array we want to know the dimensions of.
   */
  explicit Dimensions(const mxArray *ptr);

  bool are_1d() const { return are_nd(1); }//!< Is this 1D?
  bool are_2d() const { return are_nd(2); }//!< Is this 2D?
};
