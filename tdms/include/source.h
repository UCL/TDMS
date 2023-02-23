/**
 * @file source.h
 */
#pragma once

#include <complex>
#include <string>

#include "cell_coordinate.h"
#include "mat_io.h"

/**
 * @brief Typedef for providing indices to the Source class.
 *
 * Assume the notation as in the docstring of the Source class. Recall that data
 * in a Source instance is indexed via
 * {real,imag}[cell_c][cell_b][split_field_ID].
 *
 * SourceIndex is a container for indexing the Source.{real,imag} objects in an
 * efficient and notationally-consistent manner.
 */
struct SourceIndex {
  int split_field_ID = 0,//!< Index of the split field
          cell_b = 0,    //!< Index in the minor (B-)axis
          cell_c = 0;    //!< Index in the major (C-)axis
};

/**
 * @brief The Source class stores values of the Source field across a particular
 * plane.
 *
 * Let A (= {i,j,k}) be the axial direction that the plane the given Source
 * instance is storing data for. Let B, C be the remaining axial directions,
 * with C being the axial direction with the slower-varying index. That is: A =
 * i : B = j : C = k, A = j : B = i : C = k, A = k : B = i : C = j.
 *
 * The Source data (real and imag) is indexed by 3 indices, accessed via
 * {real,imag}[cell_c][cell_b][split_field_ID].
 *
 * split_field_ID ranges between 0-7 inclusive.
 * TODO: Indices <-> sources need to be deduced from input file generator
 * functions.
 *
 * Let cell_a be the A-index of the plane that the instance is storing data
 * on. Then (cell_A, cell_B, cell_C) is the Yee-cell index whose (source) data
 * we are accessing with this call.
 */
class Source {
private:
  /*! Flags if the array is empty to avoid pointer preservation */
  bool no_data_stored = true;

public:
  double ***real = nullptr;//!< Real data for the source term
  double ***imag = nullptr;//!< Imag data for the source term

  Source(const mxArray *ptr, int dim1, int dim2, const std::string &name);

  /** @brief Check if the source term is empty (true) or not (false) */
  bool is_empty() { return no_data_stored; }

  std::complex<double> operator[](SourceIndex index) {
    return std::complex<double>(
            real[index.cell_c][index.cell_b][index.split_field_ID],
            imag[index.cell_c][index.cell_b][index.split_field_ID]);
  }
};
