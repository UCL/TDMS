#include "hdf5_io.h"

#include <exception>
#include <iostream>

#include <H5Cpp.h>
#include <spdlog/spdlog.h>

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

std::vector<hsize_t> HDF5Base::shape_of(const std::string &dataname) const {

  // get the dataset and dataspace (contains dimensionality info)
  H5::DataSet dataset = file_->openDataSet(dataname);
  H5::DataSpace dataspace = dataset.getSpace();

  // need the rank in order to declare the vector size
  int rank = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dimensions(rank);
  dataspace.getSimpleExtentDims(dimensions.data(), nullptr);

  // vector is the size in each dimension i, j(, k)
  return dimensions;
}


bool HDF5Base::is_ok() const { return true; }
