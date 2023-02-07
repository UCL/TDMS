#include "hdf5_io.h"

#include <iostream>

#include <H5Cpp.h>

/******************************************************************************
 * HDF5Reader
 */
template<typename T>
void HDF5Reader::read(const std::string &dataset_name, T *data) const {

  // get the dataset and dataspace (contains dimensionality info)
  H5::DataSet dataset = file_->openDataSet(dataset_name);
  H5::DataSpace dataspace = dataset.getSpace();

  // need to get the number of matrix dimensions (rank) so that we can
  // dynamically allocate `dimensions`
  int rank = dataspace.getSimpleExtentNdims();
  hsize_t *dimensions = new hsize_t[rank];
  dataspace.getSimpleExtentDims(dimensions);

  // auto dimensions = shape_of(dataset_name);
  //  TODO why do we need `dimensions` at all here?

  // now get the data type
  H5::DataType datatype = dataset.getDataType();
  dataset.read(data, datatype);

  delete[] dimensions;
}

/******************************************************************************
 * HDF5Writer
 */
void HDF5Writer::write(const std::string &dataset_name, double *data, int size,
                       hsize_t *dimensions) {
  // 1D array
  hsize_t rank = 1;
  (void) size;// TODO what?

  // declare a dataspace
  H5::DataSpace dataspace(rank, dimensions);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

  // write the data to the dataset object in the file
  H5::DataSet dataset = file_->createDataSet(dataset_name, datatype, dataspace);
  dataset.write(data, H5::PredType::NATIVE_DOUBLE);
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

void HDF5Base::data_dump(const std::string &dataname) const {
  //
}

bool HDF5Base::is_ok() const {
  return true;
  //
}
