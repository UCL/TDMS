#pragma once

#include "hdf5_io/hdf5_base.h"

#include "arrays.h"
#include "arrays/cuboid.h"
#include "arrays/vector_typedefs.h"
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
   * @brief Read the dataset stored within a group into the buffer provided.
   * @details Can be used to read MATLAB structs by treating the struct as the
   * Group and field as the Dataset.
   * @tparam T C++ datatype to read data into.
   * @param group The Group within the file in which the dataset lives.
   * @param dataset The name of the dataset to fetch data from.
   * @param data The buffer into which to write the data.
   */
  template<typename T>
  void read_dataset_in_group(const std::string &group,
                             const std::string &dataset, T *data) const {
    spdlog::debug("Reading {} from file: {}", group, filename_);

    // Structs are saved as groups, so we need to fetch the group this struct is
    // contained in
    H5::Group structure_array = file_->openGroup(group);
    // Then fetch the requested data and read it into the buffer provided
    H5::DataSet requested_field = structure_array.openDataSet(dataset);
    requested_field.read(data, requested_field.getDataType());
  }

  /**
   * @brief Read the dataset stored within a group into the buffer provided,
   * resizing the vector buffer accordingly.
   * @details Can be used to read MATLAB structs by treating the struct as the
   * Group and field as the Dataset.
   * @tparam T C++ datatype to read data into.
   * @param group The Group within the file in which the dataset lives.
   * @param dataset The name of the dataset to fetch data from.
   * @param[out] data The buffer into which to write the data.
   */
  template<typename T>
  void read_dataset_in_group(const std::string &group,
                             const std::string &dataset,
                             std::vector<T> &data) const {
    spdlog::debug("Reading {} from file: {}", group, filename_);

    // Structs are saved as groups, so we need to fetch the group this struct is
    // contained in
    H5::Group structure_array = file_->openGroup(group);
    // Then fetch the requested data and read it into the buffer provided,
    // resizing the buffer if necessary
    H5::DataSet requested_field = structure_array.openDataSet(dataset);
    H5Dimension field_size(requested_field);
    int number_of_elements = field_size.number_of_elements();
    data.resize(number_of_elements);
    requested_field.read(data.data(), requested_field.getDataType());
  }

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
   * @brief Reads a 2D-dataset into a Matrix object.
   *
   * @tparam T C++ datatype of the Matrix object
   * @param dataset_name Name of the dataset to read data from
   * @param data_location Matrix object buffer
   */
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

  /**
   * @brief Read an InterfaceComponent into the buffer provided.
   *
   * @param[in] plane The plane {I,J,K}{0,1} to read from the file.
   * @param[out] ic InterfaceComponent reference to populate/overwrite.
   */
  void read(const std::string &plane, InterfaceComponent *ic) const;
  /**
   * @brief Read an InterfaceComponent from the file.
   *
   * @param plane The plane {I,J,K}{0,1} to read from the file.
   * @return InterfaceComponent corresponding to the requested plane.
   */
  InterfaceComponent read(const std::string &plane) const {
    InterfaceComponent ic;
    read(plane, &ic);
    return ic;
  }

  /**
   * @brief Read FrequencyVectors into the buffer provided.
   *
   * @param[out] f_vec FrequencyVectors reference to populate/overwrite.
   */
  void read(FrequencyVectors *f_vec) const;
  /**
   * @brief Read FrequencyVectors from the file.
   *
   * @return FrequencyVectors object containing the data from the input file.
   */
  FrequencyVectors read() const {
    FrequencyVectors f_vec;
    read(&f_vec);
    return f_vec;
  }

  /**
   * @brief
   *
   * Cuboid is just the phasorsurface array which is it's own dataset, but needs
   * to be offset by -1 b/c indexing things
   * @param cube
   */
  void read(Cuboid *cube) const;

  /**
   * @brief Read data from the file into a DispersiveMultiLayerObject
   * @details Data is read from the "dispersive_aux" group. If this group does
   * not exist, no data is written but no exception is thrown.
   *
   * If the group does exist, the alpha, beta, gamma, kappa, and sigma members
   * are populated with the corresponding data entries.
   * @param dml DispersiveMultiLayer object into which to write data.
   */
  void read(DispersiveMultiLayer *dml) const;

  void read(FrequencyExtractVector &fev, double omega_an,
            const std::string &dataset_name = "f_ex_vec") const;
};
