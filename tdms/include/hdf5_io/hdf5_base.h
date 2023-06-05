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
#include "hdf5_io/hdf5_dimension.h"

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
   * @brief No default constructor: it never makes sense to have a HDF5Base
   *        without a `file_` which means we need a `filename` passed to the
   *        constructor.
   */
  HDF5Base() = delete;

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
   * with `dataname`.
   * @param dataname The name of the data table.
   * @return H5Dimension The dimensions of the data.
   */
  H5Dimension shape_of(const std::string &dataname) const;
  /**
   * @brief Return shape/dimensionality information about array data stored
   * within a group.
   *
   * @param group_name The name of the HDF5 Group in which the data array is
   * stored.
   * @param dataname The name of the data array to check dimensions of.
   * @return The dimensions of the data.
   */
  H5Dimension shape_of(const std::string &group_name,
                       const std::string &dataname) const;

  /**
   * @brief Returns true if the object under /object_path is flagged as
   *        MATLAB_empty.
   *
   * @details Naturally, MATLAB does not save empty arrays or structs as objects
   *          with no elements or size, instead it saves them as 2-by-1 arrays
   *          with 0's populating the data. A simple comparison against the
   *          number of elements or members does not provide the correct
   *          information when attempting to determine if empty input has been
   *          passed.
   *
   *          MATLAB _does_, however. attach an attribute to any object that it
   *          marks as empty: MATLAB_empty. This is a uint8 set to 1 if the
   *          object is an empty array, so use that to distinguish.
   *
   * @param object_path Path to the object under file root; /object_path.
   * @return true The object is flagged as being the empty MATLAB array.
   * @return false Otherwise.
   */
  bool flagged_MATLAB_empty(const std::string &object_path) const;
};
