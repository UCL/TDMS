/**
 * @file hdf5_io.cpp
 * @authors Sam Cunliffe, William Graham
 * @brief Common HDF5 I/O methods abstracted to the base class.
 */
#include "hdf5_io/hdf5_base.h"

#include <iostream>

#include <H5Cpp.h>
#include <H5public.h>
#include <spdlog/spdlog.h>

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

H5Dimension HDF5Base::shape_of(const std::string &dataname) const {
  // Get the dataspace (contains dimensionality info)
  H5::DataSpace dataspace = file_->openDataSet(dataname).getSpace();
  return H5Dimension(dataspace);
}

H5Dimension HDF5Base::shape_of(const std::string &group_name,
                               const std::string &dataname) const {
  // Open the group that contains the dataset
  H5::Group group = file_->openGroup(group_name);
  // Get the DataSpace for the DataSet within the group
  H5::DataSpace dataspace = group.openDataSet(dataname).getSpace();
  return H5Dimension(dataspace);
}

bool HDF5Base::is_ok() const {
  // TODO: check for file health might be unnessicary given we've constructed
  // the object.
  return file_->isHdf5(filename_);
  // return file_->isAccessible(filename_) && file_->isHdf5(filename_);
}
