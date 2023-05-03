#include "hdf5_io/hdf5_writer.h"

using namespace std;

void HDF5Writer::write(const string &dataset_name, double *data, int size,
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
