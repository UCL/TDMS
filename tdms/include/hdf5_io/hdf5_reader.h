#pragma once

#include "hdf5_io/hdf5_base.h"

#include "arrays.h"
#include "interface.h"

/**
 * @brief Class wrapper of the reading of HDF5 format files.
 * @details Opens files in readonly and retrieves the datasets (in our case
 *          **double, but can be anything in general).
 */
class HDF5Reader : public HDF5Base {

public:
  /**
   * @brief Construct a new HDF5Reader for a named file.
   * @param filename The name of the file.
   * @throws H5::FileIException if the file can't be created.
   */
  HDF5Reader(const std::string &filename)
      : HDF5Base(filename, H5F_ACC_RDONLY) {}

  /**
   * @brief Reads a named dataset from the HDF5 file.
   * @param dataname The name of the datset to be read.
   * @param data A pointer to an array of correct size.
   */
  // template <typename T>
  // void read(const std::string &dataname, T *data) const;
  template<typename T>
  void read(const std::string &dataset_name, T *data) const {
    spdlog::debug("Reading {} from file: {}", dataset_name, filename_);

    // get the dataset and dataspace
    H5::DataSet dataset = file_->openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    // now get the data type
    dataset.read(data, dataset.getDataType());
    spdlog::trace("Read successful.");
  }

  template<typename T>
  void read_field_from_struct(const std::string &struct_name,
                              const std::string &field_name, T *data) const {
    spdlog::debug("Reading {} from file: {}", struct_name, filename_);

    // Structs are saved as groups, so we need to fetch the group this struct is
    // contained in
    H5::Group structure_array = file_->openGroup(struct_name);
    // Then fetch the requested data and read it into the buffer provided
    H5::DataSet requested_field = structure_array.openDataSet(field_name);
    requested_field.read(data, requested_field.getDataType());
  }

  template<typename T>
  void read(const std::string &dataset_name, Matrix<T> &data_location) const {
    spdlog::debug("Reading {} from file: {}", dataset_name, filename_);

    std::vector<hsize_t> dimensions = shape_of(dataset_name);
    if (dimensions.size() != 2) {
      throw std::runtime_error(
              "Cannot read " + dataset_name + " into a 2D matrix, it has " +
              std::to_string(dimensions.size()) + " dimensions");
    }
    int n_rows = dimensions[0];
    int n_cols = dimensions[1];

    SPDLOG_DEBUG("n_rows = {}; n_cols = {}", n_rows, n_cols);
    T *buff = (T *) malloc(n_rows * n_cols * sizeof(T));
    read(dataset_name, buff);

    data_location.allocate(n_rows, n_cols);
    for (unsigned int i = 0; i < n_rows; i++) {
      for (unsigned int j = 0; j < n_cols; j++) {
        data_location[i][j] = buff[i * n_cols + j];
      }
    }
    return;
  }

  void read(const std::string &plane, InterfaceComponent *ic) const;
  InterfaceComponent read(const std::string &plane) const {
    InterfaceComponent ic;
    read(plane, &ic);
    return ic;
  }

  void read(FrequencyVectors *f_vec) const;
  FrequencyVectors read() const {
    FrequencyVectors f_vec;
    read(&f_vec);
    return f_vec;
  }
};
