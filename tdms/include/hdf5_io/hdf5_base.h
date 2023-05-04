/**
 * @file hdf5_io.h
 * @brief Helper classes for HDF5 file I/O.
 * @details The main classes are `HDF5Reader` and `HDF5Writer` with the methods
 * `HDF5Reader::read` and `HDF5Writer::write` respectively.
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "cell_coordinate.h"

/**
 * @brief Convert from a vector of HDF5's hsize_t back to our struct of ints.
 * @note Local scope utility function as only this code needs to interact with
 * the HDF5 H5Cpp library.
 *
 * @param dimensions a 1, 2, or 3 element vector of dimensions.
 * @return ijk The dimensions in a struct.
 */
ijk to_ijk(const std::vector<hsize_t> dimensions);

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
   * @return IJKDimensions The dimensions of the data.
   */
  std::vector<hsize_t> shape_of(const std::string &dataname) const;
  /**
   * @brief Return shape/dimensionality information about the array data stored
   * within a group.
   *
   * @param group_name The name of the HDF5 Group in which the data array is
   * stored.
   * @param dataname The name of the data array to check dimensions of.
   * @return std::vector<hsize_t> The dimensions of the data.
   */
  std::vector<hsize_t> shape_of(const std::string &group_name,
                                const std::string &dataname) const;

  /**
   * @brief Checks the file is a valid HDF5 file, and everything is OK.
   * TODO: Can perhaps remove.
   *
   * @return true If all is well.
   * @return false Otherwise.
   */
  bool is_ok() const;
};
