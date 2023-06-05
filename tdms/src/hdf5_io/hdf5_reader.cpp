#include "hdf5_io/hdf5_reader.h"
#include "hdf5_io/hdf5_dimension.h"

#include <algorithm>
#include <stdexcept>

#include <spdlog/spdlog.h>

using namespace std;

void HDF5Reader::read(const string &plane, InterfaceComponent *ic) const {
  // Read the InterfaceComponent in as a 2-element double array
  double read_buffer[2];
  read_dataset_in_group("interface", plane, read_buffer);
  // The index that is read in should have 1 subtracted from it, to account for
  // MATLAB indexing
  ic->index = max((int) read_buffer[0] - 1, 0);
  // The apply flag should be cast from the double that is read in
  ic->apply = (bool) read_buffer[1];
}

void HDF5Reader::read(FrequencyVectors *f_vec) const {
  // Allocate memory in f_vec
  H5Dimension x_dims = shape_of("f_vec/fx_vec");
  H5Dimension y_dims = shape_of("f_vec/fy_vec");
  // Check that we have one dimensional arrays
  if (!x_dims.is_1D() || !y_dims.is_1D()) {
    throw runtime_error("f_vec members are not 1D arrays!");
  }
  // Allocate memory - resize() must be used over reserve() to ensure enough
  // buffer when we read in, _and_ that size() correctly returns the number of
  // elements in the buffer.
  f_vec->x.resize(x_dims.max_dim());
  f_vec->y.resize(y_dims.max_dim());
  // Now read the data into the vectors
  read_dataset_in_group("f_vec", "fx_vec", f_vec->x.data());
  read_dataset_in_group("f_vec", "fy_vec", f_vec->y.data());
}

void HDF5Reader::read(Cuboid *cube) const {
  string cuboid_dataset = "phasorsurface";
  // Check that we are reading in a 1D array with 6 elements
  H5Dimension cuboid_dims(file_->openDataSet(cuboid_dataset).getSpace());
  if (!(cuboid_dims.is_1D() && cuboid_dims.max_dim() == 6)) {
    throw runtime_error(
            "Error: phasorsurface is not a 1D vector of 6 elements");
  }
  // Read buffer then adjust for indexing offset between MATLAB and C++
  // NOTE: Buffer is saved as doubles in .mat file, but we want to read as
  // integers here.
  double intermediate_buffer[6];
  read(cuboid_dataset, intermediate_buffer);
  for (int i = 0; i < 6; i++) {
    cube->array[i] = (int) intermediate_buffer[i] - 1;
  }
}

void HDF5Reader::read(DispersiveMultiLayer *dml) const {
  string group_name = "dispersive_aux";
  // Deal with the case of an empty input
  if (!file_->nameExists(group_name)) {
    spdlog::info(group_name + " is not a group: assuming empty input");
    return;
  } else {
    // This is a group - it should have 9 members and we can quickly check this
    H5::Group dispersive_aux = file_->openGroup(group_name);
    if (dispersive_aux.getNumObjs() != 9) {
      throw runtime_error("dispersive_aux does not have exactly 9 members!");
    }
  }
  // Assuming non-empty input, setup the data appropriately
  read_dataset_in_group(group_name, "alpha", dml->alpha);
  read_dataset_in_group(group_name, "beta", dml->beta);
  read_dataset_in_group(group_name, "gamma", dml->gamma);
  read_dataset_in_group(group_name, "kappa_x", dml->kappa.x);
  read_dataset_in_group(group_name, "kappa_y", dml->kappa.y);
  read_dataset_in_group(group_name, "kappa_z", dml->kappa.z);
  read_dataset_in_group(group_name, "sigma_x", dml->sigma.x);
  read_dataset_in_group(group_name, "sigma_y", dml->sigma.y);
  read_dataset_in_group(group_name, "sigma_z", dml->sigma.z);
}
