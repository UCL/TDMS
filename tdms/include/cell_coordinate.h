/**
 * @file cell_coordinate.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Class declaration for Yee cell coordinates
 */
#pragma once

class CellCoordinate {

private:
  // the indices of the Yee cell. Private so that they cannot be edited unintentionally
  int cell_i = 0, cell_j = 0, cell_k = 0;

public:
  /**
   * @brief Construct a new Cell Coordinate object
   *
   * @param i,j,k The cell index to assign
   */
  CellCoordinate(int i = 0, int j = 0, int k = 0) {
    cell_i = i;
    cell_j = j;
    cell_k = k;
  };
  /**
   * @brief Construct a new Cell Coordinate object, passing in a contiguous block of memory
   *
   * @param indices The indices of the cell index to assign, in the order i,j,k
   * @param buffer_start First index to read the indices from
   */
  CellCoordinate(int *indices, int buffer_start = 0) {
    cell_i = indices[buffer_start];
    cell_j = indices[buffer_start + 1];
    cell_k = indices[buffer_start + 2];
  };

  int i() { return cell_i; };   //< fetch i
  int j() { return cell_j; };   //< fetch j
  int k() { return cell_k; };   //< fetch k
};
