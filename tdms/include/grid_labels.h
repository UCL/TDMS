/**
 * @file grid_labels.h
 * @brief Class to hold the labels of a Yee cell.
 */
#pragma once

#include "mat_io.h"

/**
 * Grid labels hold the cartesian labels of Yee cell, in the x, y and z
 * directions
 */
class GridLabels {
public:
  double *x = nullptr;//< Start of the labels in the x direction
  double *y = nullptr;//< Start of the labels in the y direction
  double *z = nullptr;//< Start of the labels in the z direction

  GridLabels() = default;

  explicit GridLabels(const mxArray *ptr);

  /**
   * @brief Set values by copying from another GridLabels object
   *
   * @param other_labels The GridLabels object to copy values from
   * @param i_l,j_l,k_l The first item to copy from the x, y, z attributes
   * (respectively)
   * @param i_u,j_u,k_u The final (inclusive) item to copy from the x, y, z
   * attributes (respectively)
   */
  void initialise_from(const GridLabels &labels_to_copy_from, int i_l, int i_u,
                       int j_l, int j_u, int k_l, int k_u);
};
