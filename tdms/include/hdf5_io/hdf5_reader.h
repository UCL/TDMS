#pragma once

#include "hdf5_io/hdf5_base.h"
#include "hdf5_io/hdf5_dimension.h"

#include "arrays.h"
#include "arrays/cuboid.h"
#include "arrays/tdms_matrix.h"
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

  /**
   * @brief Reads the data from the dataset at the location (under file root)
   * provided.
   * @tparam T Datatype to read.
   * @param path_to_dataset Path under file root of the dataset to read from.
   * @param buffer Vector to read data into.
   */
  template<typename T>
  void read(const std::string &path_to_dataset, std::vector<T> &buffer) const {
    spdlog::debug("Reading {} from file", path_to_dataset, filename_);

    if (path_exists(path_to_dataset, H5I_DATASET)) {
      H5::DataSet dataset = file_->openDataSet(path_to_dataset);
      H5Dimension dims = shape_of(path_to_dataset);
      buffer.resize(dims.number_of_elements());
      dataset.read(buffer.data(), dataset.getDataType());
    } else {
      // Path does not exist, throw error
      throw std::runtime_error(path_to_dataset +
                               " does not point to a dataset");
    }
  }

  /**
   * @brief Reads a 2D-dataset into a Matrix object.
   *
   * @tparam T C++ datatype of the Matrix object
   * @param dataset_name Name of the dataset to read data from
   * @param data_location Matrix object buffer
   */
  template<typename T>
  void read(const std::string &dataset_name, Matrix<T> &data_location) const {
    spdlog::debug("Reading {} from file: {}", dataset_name, filename_);

    H5Dimension dimensions = shape_of(dataset_name);
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
        data_location(i, j) = buff[i * n_cols + j];
      }
    }
    return;
  }

  /**
   * @brief Read an InterfaceComponent into the buffer provided.
   *
   * @param[in] plane The plane {I,J,K}{0,1} to read from the file.
   * @param[out] ic InterfaceComponent reference to populate/overwrite.
   */
  void read(const std::string &plane, InterfaceComponent &ic) const;

  /**
   * @brief Read FrequencyVectors into the buffer provided.
   *
   * @param[out] f_vec FrequencyVectors reference to populate/overwrite.
   */
  void read(FrequencyVectors &f_vec) const;

  /**
   * @brief
   *
   * Cuboid is just the phasorsurface array which is it's own dataset, but needs
   * to be offset by -1 b/c indexing things
   * @param cube
   */
  void read(Cuboid &cube) const;

  /**
   * @brief Read data from the file into a DispersiveMultiLayerObject
   * @details Data is read from the "dispersive_aux" group. If this group does
   * not exist, no data is written but no exception is thrown.
   *
   * If the group does exist, the alpha, beta, gamma, kappa, and sigma members
   * are populated with the corresponding data entries.
   * @param dml DispersiveMultiLayer object into which to write data.
   */
  void read(DispersiveMultiLayer &dml) const;
};
