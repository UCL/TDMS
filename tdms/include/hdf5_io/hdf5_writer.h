#pragma once

#include "hdf5_io/hdf5_base.h"

#include "arrays/tdms_matrix.h"

class HDF5Writer : public HDF5Base {

public:
  /**
   * @brief Construct a new HDF5Writer, creates a file.
   * @param filename The name of the file to be created.
   * @throws H5::FileIException if the file can't be created.
   */
  HDF5Writer(const std::string &filename) : HDF5Base(filename, H5F_ACC_TRUNC) {}

  /**
   * @brief Write `data` to the file with `dataname`.
   *
   * @param dataname The name of the data table.
   * @param data The data itself.
   * @param size The size of the data array.
   * @param dimensions The number of dimensions of the array.
   */
  void write(const std::string &dataname, const double *data, int size,
             hsize_t *dimensions);

  /**
   * @brief Write `data` to the file with `dataname`.
   *
   * @param dataname The name of the data table.
   * @param data The data itself.
   * @param size The size of the data array.
   * @param dimensions The number of dimensions of the array.
   */
  template<typename T>
  void write(const std::string &dataname, const Matrix<T> &data) {
    hsize_t dimensions[2];
    dimensions[0] = data.get_n_rows();
    dimensions[1] = data.get_n_cols();
    write(dataname, data.data_.data(), 2, dimensions);
  }
};
