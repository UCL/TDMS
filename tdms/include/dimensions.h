/**
 * @file dimensions.h
 * @brief Defines a class that serves as an explicit converter between MATLAB
 * pointers and array dimensions
 */
#pragma once

#include "mat_io.h"

class Dimensions {
private:
  int i = 0;//!< Extent of dimension 1
  int j = 0;//!< Extent of dimension 2
  int k = 0;//!< Extent of dimension 3

  bool are_nd(int n) const { return (bool(i) + bool(j) + bool(k)) == n; }

public:
  int operator[](int value) const;

  explicit Dimensions(const mxArray *ptr);

  bool are_1d() const { return are_nd(1); }
  bool are_2d() const { return are_nd(2); }
};
