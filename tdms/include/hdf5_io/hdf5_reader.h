#pragma once

#include <stdexcept>
#include <string>

#include "hdf5_io/hdf5_base.h"
#include "hdf5_io/hdf5_dimension.h"

#include "arrays.h"
#include "arrays/cuboid.h"
#include "arrays/dispersive_multilayer.h"
#include "arrays/material_collections.h"
#include "arrays/tdms_matrix.h"
#include "arrays/vector_typedefs.h"
#include "arrays/xyz_vector.h"
#include "interface.h"

/**
 * @brief Class wrapper of the reading of HDF5 format files.
 * @details Opens files in readonly and retrieves the datasets (in our
 *          case, **double, but can be anything in general).
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

  /**
   * @brief Read data from a group with {x,y,z} members containing 1D data.
   * @details There is no obligation for the three members to be read from to be
   * of the same size, however they are expected to be doubles.
   *
   * The datasets themselves do not have to literally be called x, y, and z, but
   * they must all be named with a common prefix that terminates in one of these
   * three characters. This prefix must be specified by passing the name_prefix
   * argument (an empty string is an acceptable value to pass).
   *
   * For example,
   * group
   * - data_x
   * - data_y
   * - data_z
   * would be acceptable, however
   * group
   * - x_data
   * - y_data
   * - z_data
   * and
   * group
   * - data_x
   * - data_for_y
   * - data_z
   * would be unacceptable.
   * @param group_name Name of the group within the file to read from.
   * @param name_prefix A prefix to be attached to the {x,y,z} names of the
   * datasets.
   * @param v The XYZVector to read the data into.
   */
  void read(const std::string &group_name, const std::string &name_prefix,
            XYZVector &v) const;

  /**
   * @brief Read C-constants (algebraic terms) from the file.
   * @details See CMaterial docstring for more information about the constants
   * to be read.
   * @param c_material Structure to read information into.
   * @param group_name The name of the group in the file to read the data from.
   * Defaults to Cmaterial (the expected location) if not provided.
   */
  void read(CMaterial &c_material,
            const std::string &group_name = "Cmaterial") const;

  /** @overload void read(CMaterial &c_material, const std::string &group_name =
   * "Cmaterial")
   */
  void read(CCollection &c_collection,
            const std::string &dataset_name = "C") const;

  /**
   * @brief Read D-constants (algebraic terms) from the file.
   * @details See DBase docstring for more information about the constants to be
   * read.
   * @tparam is_material Either true or false, indicating either a DMaterial or
   * DCollection is to be read. These structures are identical in terms of
   * memory, however the default read location in the file differs between the
   * two.
   * @param d_base Either a DMaterial or DCollection, the structure to read the
   * data into.
   * @param group_name The name of the group in the file to read the data from.
   * Uses the default location if not provided.
   */
  template<bool is_material>
  void read(DBase<is_material> &d_base,
            const std::string group_name = "") const {
    // Read from non-standard location in file if passed a second argument
    const std::string group_to_read_from =
            group_name.empty() ? d_base.input_field : group_name;

    // We should expect a group with 6 members
    H5::Group group = file_->openGroup(group_to_read_from);
    int n_members = group.getNumObjs();
    if (n_members != 6) {
      throw std::runtime_error("D should have 6 members, but " +
                               std::to_string(n_members) + " were found");
    }
    group.close();

    read(group_to_read_from, "Da", d_base.a);
    read(group_to_read_from, "Db", d_base.b);
  };

  /**
   * @brief Read data from the file into a FrequencyExtractVector.
   * @details If the data to be read is one-dimensional and non-empty, the
   * values are simply read as as usual. If the dataset that is pointed to is
   * _empty_, then the FrequencyExtractVector is set to contain a single,
   * default value of \f$\omega_{an}/2\pi\f$.
   * @param fev The FrequencyExtractVector to place the data into.
   * @param omega_an Angular frequency.
   * @param dataset_name Dataset within the file to read from.
   */
  void read(FrequencyExtractVector &fev, double omega_an,
            const std::string &dataset_name = "f_ex_vec") const;
};
