#include "hdf5_io/hdf5_reader.h"
#include "hdf5_io/hdf5_dimension.h"

#include <algorithm>
#include <stdexcept>

using namespace std;

void HDF5Reader::read(const string &plane, InterfaceComponent *ic) const {
  // Read the InterfaceComponent in as a 2-element double array
  double read_buffer[2];
  read_field_from_struct("interface", plane, read_buffer);
  // The index that is read in should have 1 subtracted from it, to account for
  // MATLAB indexing
  ic->index = max((int) read_buffer[0] - 1, 0);
  // The apply flag should be cast from the double that is read in
  ic->apply = (bool) read_buffer[1];
}

void HDF5Reader::read(FrequencyVectors *f_vec) const {
  // Allocate memory in f_vec
  H5Dimension x_dims = shape_of("f_vec", "fx_vec");
  H5Dimension y_dims = shape_of("f_vec", "fy_vec");
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
  read_field_from_struct("f_vec", "fx_vec", f_vec->x.data());
  read_field_from_struct("f_vec", "fy_vec", f_vec->y.data());
}
