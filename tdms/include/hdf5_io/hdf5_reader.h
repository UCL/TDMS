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
   * @brief Read data from the file, expecting the given data type to be present
   * in the file.
   * @details To read non-standard datatypes from files, we must define a
   * corresponding structure array that matches the data-type, and a
   * H5::CompType from this which can be used to read the data in via the HDF5
   * API.
   *
   * The compound type should be created via
   * H5::CompType new_comp_type(sizeof(T));
   * and then have the data-type members added via the
   * new_comp_type.insertMember() method.
   * @tparam T Datatype represented in HDF5 by the parameter compound_type.
   * @param path_to_dataset Path from the file root to the dataset to read.
   * @param[out] buffer Buffer into which to read the data. Will be
   * automatically resized.
   * @param compound_type H5::CompType that the data is stored as, in the
   * dataset.
   */
  template<typename T>
  void read(const std::string &path_to_dataset, std::vector<T> &buffer,
            const H5::CompType &compound_type) const {
    spdlog::debug("Reading COMPOUND DATA from {} in {}", path_to_dataset,
                  filename_);

    if (path_exists(path_to_dataset, H5I_DATASET)) {
      H5::DataSet dataset = file_->openDataSet(path_to_dataset);
      H5Dimension dims = shape_of(path_to_dataset);
      buffer.resize(dims.number_of_elements());
      dataset.read(buffer.data(), compound_type);
    } else {
      // Path does not exist, throw error
      throw std::runtime_error(path_to_dataset + "does not point to a dataset");
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

    data_location.allocate(dimensions[0], dimensions[1]);
    read(dataset_name, data_location.data_.data());
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
