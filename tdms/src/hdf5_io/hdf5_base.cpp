/**
 * @file hdf5_io.cpp
 * @authors Sam Cunliffe, William Graham
 * @brief Common HDF5 I/O methods abstracted to the base class.
 */
#include "hdf5_io/hdf5_base.h"

#include <iostream>
#include <stdexcept>

#include <H5Cpp.h>
#include <H5public.h>
#include <spdlog/spdlog.h>

bool HDF5Base::path_exists(const std::string &path_under_root) const {
  // Can't check anything if there's no file
  if (file_ == nullptr) { throw std::runtime_error("No file opened"); }

  // Attempt to lookup the path
  return file_->exists(path_under_root);
}

bool HDF5Base::path_exists(const std::string &path_under_root,
                           H5I_type_t *points_to) const {
  // Attempt to lookup the path
  bool return_value = path_exists(path_under_root);

  // Set object type if the path existed
  if (return_value) {
    hid_t object_reference = file_->getObjId(path_under_root);
    *points_to = file_->getHDFObjType(object_reference);
  }

  return return_value;
}

bool HDF5Base::path_exists(const std::string &path_under_root,
                           const H5I_type_t &object_type,
                           const bool error_on_false) const {
  // Attempt to lookup the path
  bool path_is_valid = path_exists(path_under_root);
  bool object_is_correct_type = false;
  // If the path is valid, check the object it points to is the correct type
  if (path_is_valid) {
    hid_t object_reference = file_->getObjId(path_under_root);
    object_is_correct_type =
            file_->getHDFObjType(object_reference) == object_type;
  }
  // Return result, or throw error if running strictly
  if (path_is_valid && object_is_correct_type) {
    return true;
  } else if (error_on_false) {
    throw std::runtime_error(path_under_root +
                             "does not point to an object of the correct type");
  }
  return false;
}

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

H5Dimension HDF5Base::shape_of(const std::string &dataset_path) const {
  if (path_exists(dataset_path, H5I_type_t::H5I_DATASET)) {
    // Get the dataspace (contains dimensionality info)
    return H5Dimension(file_->openDataSet(dataset_path).getSpace());
  } else {
    throw std::runtime_error(dataset_path + " does not point to a dataset");
  }
}

bool HDF5Base::is_ok() const {
  // TODO: check for file health might be unnecessary given we've constructed
  // the object.
  return file_->isHdf5(filename_);
  // return file_->isAccessible(filename_) && file_->isHdf5(filename_);
}

bool HDF5Base::flagged_MATLAB_empty(const std::string &object_path) const {
  H5::Attribute empty_attribute;// will point to the MATLAB_empty attribute

  if (path_exists(object_path, H5I_GROUP)) {
    // Dealing with a group
    H5::Group object = file_->openGroup(object_path);
    if (object.attrExists("MATLAB_empty")) {
      empty_attribute = object.openAttribute("MATLAB_empty");
    } else {
      // Object is not flagged as being empty
      return false;
    }
    object.close();
  } else if (path_exists(object_path, H5I_DATASET)) {
    // Dealing with a dataset
    H5::DataSet object = file_->openDataSet(object_path);
    if (object.attrExists("MATLAB_empty")) {
      empty_attribute = object.openAttribute("MATLAB_empty");
    } else {
      // Object is not flagged as being empty
      return false;
    }
    object.close();
  } else {
    // No other objects should be the result of MATLAB saving an empty object,
    // so throw error
    throw std::runtime_error(object_path + " is not a Group or a DataSet");
  }

  // Having extracted the dataset, attempt to read the attribute value
  uint8_t empty_bool[1];
  empty_attribute.read(H5::PredType::NATIVE_UINT8, empty_bool);
  empty_attribute.close();

  // And finally return whether this is flagged as true
  return empty_bool[0] == 1;
}
