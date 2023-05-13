/**
 * @file hdf5_io.h
 * @brief Helper classes for HDF5 file I/O.
 * @details The main classes are `HDF5Reader` and `HDF5Writer` with the methods
 * `HDF5Reader::read` and `HDF5Writer::write` respectively.
 */
#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "hdf5_io/hdf5_dimension.h"

/**
 * @brief Convert from C++ datatype T to a H5::PredType for HDF5 write out.
 *
 * TODO: Doesn't support complex doubles. Might need a compount datatype.
 * Investigate how to read complex doubles first.
 *
 * @tparam T the native C++ type.
 * @return H5::PredType The HDF5 type to save.
 */
template<typename T>
H5::PredType to_hdf5_datatype() {

  if (std::is_same<T, double>()) {
    return H5::PredType::NATIVE_DOUBLE;
  } else if (std::is_same<T, int>()) {
    return H5::PredType::NATIVE_INT;
  } else if (std::is_same<T, char>()) {
    return H5::PredType::NATIVE_CHAR;
  } else if (std::is_same<T, uint16_t>()) {
    return H5::PredType::NATIVE_UINT16;
  } else if (std::is_same<T, int16_t>()) {
    return H5::PredType::NATIVE_INT16;
  } else {
    throw std::runtime_error("I don't know how to convert this datatype to "
                             "something HDF5 wants");
  }
};

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
   * with `dataname`.
   * @param dataname The name of the data table.
   * @return The dimensions of the data.
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
   * @brief Checks the file is a valid HDF5 file, and everything is OK.
   * TODO: Can perhaps remove.
   *
   * @return true If all is well.
   * @return false Otherwise.
   */
  bool is_ok() const;

  /**
   * @brief Checks if `name` exists in this file.
   *
   * @details Just wraps some ugly H5Lexists call from the C-API.
   *          TODO: does this work if given 'group/subgroup' ?
   *
   * @param name Name of the H5 object (file or group) under root to search for
   * @return true The object is present in the file, under root
   * @return false Otherwise
   */
  bool contains(const std::string &name) const {
    return H5Lexists(file_->getId(), name.c_str(), H5P_DEFAULT) > 0;
  };
};
