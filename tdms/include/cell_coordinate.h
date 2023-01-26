/**
 * @file cell_coordinate.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Container for {ijk} or {IJK} grouped variables
 */
#pragma once

#include <spdlog/spdlog.h>

#include <algorithm>


/**
 * @brief A structure for holding three values, which typically pertain to the same quantity but for each of the axial directions.
 *
 * Effectively stores the 3-vector (i,j,k). This is typically used to represent Yee cell indices, or 3D-array dimensions, or the maximum number of Yee cells in each coordinate direction, for example.
 */
struct ijk {
  // The values (by dimension) that the object contains
  int i = 0, j = 0, k = 0;

  /** @brief Return the maximum of i,j,k */
  int max() const { return std::max(std::max(i, j), k); }

  /** @brief Print the (i,j,k) values */
  void print() const { spdlog::info("ijk: ({},{},{})", i, j, k); }

  ijk &operator+=(int n) {
    this->i += n;
    this->j += n;
    this->k += n;
    return *this;
  }
};

/* Synonyms for code readability */

typedef ijk CellCoordinate;//< Index-coordinates (i,j,k) of a Yee cell
typedef ijk IJKDimensions; //< Holds array dimensions (I_tot, J_tot, K_tot)
