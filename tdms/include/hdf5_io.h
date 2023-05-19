/**
 * @file hdf5_io.h
 * @brief Helper classes for HDF5 file I/O.
 * @details The main classes are `HDF5Reader` and `HDF5Writer` with the methods
 * `HDF5Reader::read` and `HDF5Writer::write` respectively.
 */
#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <H5Cpp.h>
#include <spdlog/spdlog.h>

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
   * @throws XX if file is not found.
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

    // get the dataset and dataspace (contains dimensionality info)
    H5::DataSet dataset = file_->openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();
    spdlog::debug("Created dataspace");

    // need to get the number of matrix dimensions (rank) so that we can
    // dynamically allocate `dimensions`
    int rank = dataspace.getSimpleExtentNdims();
    spdlog::debug("Rank of dataspace: {}", rank);
    hsize_t *dimensions = new hsize_t[rank];
    dataspace.getSimpleExtentDims(dimensions);
    spdlog::debug("Got dimensions");

    // now get the data type
    H5::DataType datatype = dataset.getDataType();
    spdlog::debug("Got datatype");
    dataset.read(data, datatype);
    spdlog::debug("Read");

    delete[] dimensions;
  }
};

class HDF5Writer : public HDF5Base {

public:
  /**
   * @brief Construct a new HDF5Writer, creates a file.
   * @param filename The name of the file to be created.
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
};
