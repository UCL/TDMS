/**
 * @file hdf5_io.h
 * @brief Helper classes for HDF5 file I/O.
 * @details The main classes are `HDF5Reader` and `HDF5Writer` with the methods
 * `HDF5Reader::read` and `HDF5Writer::write` respectively.
 */
#pragma once

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <H5Cpp.h>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "fdtd_grid_initialiser.h"

/**
 * @brief The base class for HDF5 I/O.
 * @details Common functionality and wraps handling the std::unique_ptr to hold
 * the H5::File object.
 */
class HDF5Base {

protected:
  std::string filename_;             /**< The name of the file. */
  std::shared_ptr<H5::H5File> file_; /**< Pointer to the underlying H5::File. */

  /**
   * @brief Construct a new HDF5{Reader/Writer} for a named file.
   * @param filename The name of the file.
   * @param mode The H5 file access mode (RDONLY for a HDF5Reader, TRUNC for a
   * HDF5Writer.)
   * @throws H5::FileIException if the file doesn't exist or can't be created.
   */
  HDF5Base(const std::string &filename, int mode = H5F_ACC_RDONLY)
      : filename_(filename) {
    file_ = std::make_unique<H5::H5File>(filename, mode);
  }

  /**
   * @brief Destructor closes the file.
   * @details Closes file when HDF5Reader(or HDF5Writer) goes out of scope.
   * Since the file pointer is a smart pointer it is deallocated automatically.
   */
  ~HDF5Base() { file_->close(); }

public:
  /**
   * @brief Get the name of the file.
   * @return std::string the filename.
   */
  std::string get_filename() const { return filename_; }

  /**
   * @brief Get the names of all datasets (data tables) currently in the file.
   * @return std::vector<std::string> A vector of their names.
   */
  std::vector<std::string> get_datanames() const;

  /**
   * @brief Print the names of all datasets to std::out.
   */
  void ls() const;

  /**
   * @brief Return shape/dimensionality information about the array data stored
   * with `name`.
   * @param dataname The name of the data table.
   * @return std::vector<hsize_t> The dimensions of the data.
   */
  // IJKDimensions shape_of(const std::string &dataname) const;
  std::vector<hsize_t> shape_of(const std::string &dataname) const;

  /**
   * @brief Checks the file is a valid HDF5 file, and everything is OK.
   * TODO: Can perhaps remove.
   *
   * @return true If all is well.
   * @return false Otherwise.
   */
  bool is_ok() const;
};

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

  void read(const fdtdGridInitialiser &initialiser,
            const std::string &dataset_name = "fdtdgrid") const {

    // This method will take in the fdtdGridInit... object, assign the pointer
    // member, and then run the object's
    throw std::logic_error("Not yet implemented");
    return;
  }
};

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
};
