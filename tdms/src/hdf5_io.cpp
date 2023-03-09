#include "hdf5_io.h"
#include "cell_coordinate.h"

#include <exception>
#include <iostream>

#include <H5Cpp.h>
#include <H5public.h>
#include <spdlog/spdlog.h>

/**
 * @brief Convert from a vector of HDF5's hsize_t back to our struct of ints.
 * @note Local scope utility function as only this code needs to interact with
 * the HDF5 H5Cpp library.
 *
 * @param dimensions a 1, 2, or 3 element vector of dimensions.
 * @return ijk The dimensions in a struct.
 */
ijk to_ijk(const std::vector<hsize_t> dimensions) {
  unsigned int rank = dimensions.size();
  ijk out;
  if (rank > 0) out.i = (int) dimensions[0];
  if (rank > 1) out.j = (int) dimensions[1];
  if (rank > 2) out.k = (int) dimensions[2];
  if (rank > 3) spdlog::warn("Rank > 3");
  return out;
}

/******************************************************************************
 * HDF5Writer
 */
void HDF5Writer::write(const std::string &dataset_name, double *data, int size,
                       hsize_t *dimensions) {
  spdlog::debug("Writing {} to file: {}", dataset_name, filename_);

  // declare a dataspace
  H5::DataSpace dataspace(size, dimensions);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

  // write the data to the dataset object in the file
  H5::DataSet dataset = file_->createDataSet(dataset_name, datatype, dataspace);
  dataset.write(data, H5::PredType::NATIVE_DOUBLE);
  spdlog::trace("Write successful.");
}

/******************************************************************************
 * HDF5Base
 *
 * Common HDF5 I/O methods abstracted to the base class.
 */
std::vector<std::string> HDF5Base::get_datanames() const {
  std::vector<std::string> names;

  // iterate over all objects in the file
  for (unsigned int i = 0; i < file_->getNumObjs(); i++) {
    H5G_obj_t object_type = file_->getObjTypeByIdx(i);

    // if the current object is a H5::Dataset then grab its name
    if (object_type == H5G_DATASET) {
      H5std_string object_name = file_->getObjnameByIdx(i);
      names.push_back(object_name);
    }
  }
  return names;
}

void HDF5Base::ls() const {
  std::vector<std::string> names = this->get_datanames();
  for (auto name : names) std::cout << name << std::endl;
  return;
}

// std::vector<hsize_t> HDF5Base::shape_of(const std::string &dataname) const {
IJKDimensions HDF5Base::shape_of(const std::string &dataname) const {

  // get the dataset and dataspace (contains dimensionality info)
  H5::DataSet dataset = file_->openDataSet(dataname);
  H5::DataSpace dataspace = dataset.getSpace();

  // need the rank in order to declare the vector size
  int rank = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dimensions(rank);
  dataspace.getSimpleExtentDims(dimensions.data(), nullptr);
  // vector is the size in each dimension i, j(, k)
  return to_ijk(dimensions);
}

bool HDF5Base::is_ok() const {
  // TODO: check for file health might be unnessicary fiven we've constructed
  // the object.
  return file_->isAccessible(filename_) && file_->isHdf5(filename_);
}
