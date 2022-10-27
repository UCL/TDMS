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
  CellCoordinate(int i = 0, int j = 0, int k = 0) {
    cell_i = i;
    cell_j = j;
    cell_k = k;
  };
  CellCoordinate(int *indices, int buffer_start = 0) {
    cell_i = indices[buffer_start];
    cell_j = indices[buffer_start + 1];
    cell_k = indices[buffer_start + 2];
  };

  int i() { return cell_i; };   //< fetch i
  int j() { return cell_j; };   //< fetch j
  int k() { return cell_k; };   //< fetch k

  /**
   * @brief Copies the values of the cell indices into the addresses provided
   *
   * @param[out] i,j,k Address to place the value of cell_{i,j,k} into (respectively)
   */
  void unpack(int *i, int *j, int *k) {
    *i = cell_i;
    *j = cell_j;
    *k = cell_k;
  }
};
