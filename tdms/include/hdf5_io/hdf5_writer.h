#pragma once

#include <H5Lpublic.h>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "hdf5_io/hdf5_base.h"

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
  void write(const std::string &dataname, double *data, int size,
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
    int n_cols = data.get_n_cols();
    int n_rows = data.get_n_rows();
    hsize_t dimension[2] = {static_cast<hsize_t>(n_rows),
                            static_cast<hsize_t>(n_cols)};
    T *buff = (T *) malloc(n_rows * n_cols * sizeof(T));
    for (unsigned int i = 0; i < n_rows; i++) {
      for (unsigned int j = 0; j < n_cols; j++) {
        buff[i * n_cols + j] = data[i][j];
      }
    }
    write(dataname, buff, 2, dimension);
  }

  /**
   * @brief Write the dataset to a group, creating the group if it does not
   * exist.
   *
   * @details Can be used to write MATLAB structs by treating the struct as the
   *          Group and field as the Dataset.
   *
   * @tparam T      C++ datatype to write.
   * @param group   The Group to write the dataset to (will create if doesn't
   * exist).
   * @param dataset The name of the dataset to write to.
   * @param data    The buffer from which to read the data.
   */
  template<typename T>
  void write_dataset_to_group(const std::string &group,
                              const std::string &dataset,
                              std::vector<T> data) const {
    spdlog::debug("Writing {} to file: {}", group, filename_);

    // Create group if it does not exist
    H5::Group group_to_write_to;
    if (!this->contains(group)) {
      // Group does not exist, create it and open it
      group_to_write_to = file_->createGroup(group.c_str());
    } else {
      // Open the group
      // TODO how do we know this is a group and not a dataset?
      group_to_write_to = file_->openGroup(group.c_str());
    }

    // Write the data to the group - do we need to write dimensions too? Discuss
    // TODO: compile with this as unsigned long long and see if we get compiler
    // warnings.
    hsize_t dims[1] = {data.size()};
    H5::DataSpace dimension_info(1, dims);
    H5::DataSet dataset_to_write_to = group_to_write_to.createDataSet(
            dataset,
            // to_hdf5_dtype(data),
            H5::PredType::NATIVE_DOUBLE, dimension_info);

    dataset_to_write_to.write(data.data(),
                              // to_hdf5_dtype(data),
                              H5::PredType::NATIVE_DOUBLE);

    dataset_to_write_to.close();
    dimension_info.close();
    group_to_write_to.close();
    return;
  }
};

// https://support.hdfgroup.org/HDF5/doc/RM/RM_H5T.html#Datatype-Equal <-
// possible compare C++ types to HDF5 types?
