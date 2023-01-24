/**
 * @file cell_coordinate.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Container for {ijk} or {IJK} grouped variables
 */
#pragma once

#include <spdlog/spdlog.h>

struct ijk {
  // The values (by dimension) that the object contains
  int i = 0, j = 0, k = 0;

  int max() { return std::max(std::max(i, j), k); }
  void print() { spdlog::info("ijk: ({},{},{})", i, j, k); }
};

// synonyms for code readability
typedef ijk CellCoordinate;
typedef ijk IJKDims;
