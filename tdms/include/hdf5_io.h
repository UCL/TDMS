/**
 * @file hdf5_io.h
 * @brief File I/O using HDF5.
 */
#pragma once
#include <string>
#include <vector>

#include <H5Cpp.h>

/**
 * @brief Mode to write the file: R, RW, O.
 */
enum HDF5FileMode {
  READONLY = 0,
  READWRITE,
  OVERWRITE
};

/**
 * @brief A simple class to wrap the HDF5 I/O.
 */
class HDF5File {

public:
  /** Default constructor: default filename and mode. */
  HDF5File() : filename_("tdms.hdf5"), mode_(READONLY) { _open(); }

  /** Normal constructor. */
  HDF5File(const std::string &filename, HDF5FileMode mode) : filename_(filename), mode_(mode) {
    _open();
  }

  /** Destructor, deletes pointers. */
  ~HDF5File();

  /** Writes array to the file. */
  void write();

  /** reads data from the file. */
  void read();

  /**
   * @brief Check file health.
   *
   * @param print_debug Optionally this function can print to the debug log.
   * @return true If the file is OK.
   * @return false If the file is not OK.
   */
  bool isOK(bool print_debug=false);

private:
  /** Common to both constructors: open/create the file and set file_ to point to it. */
  void _open();
  std::string filename_;                /**< The file name. */
  HDF5FileMode mode_;                   /**< The I/O mode: default is to create non-existing. */
  H5::H5File *file_ = nullptr;          /**< The H5 file itself. */
  std::vector<H5::DataSet *> datasets_; /**< All datasets to be written to the file. */
};

/**
 * @brief Example of HDF5 I/O... to be deleted.
 * This function should not make it to a PR.
 */
void example_hdf5();
