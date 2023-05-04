#include "hdf5_io/hdf5_reader.h"

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
  vector<hsize_t> x_dims = shape_of("f_vec", "fx_vec");
  vector<hsize_t> y_dims = shape_of("f_vec", "fy_vec");
  // Check that we have ONE DIMENSIONAL ARRAY! FOR NOW LETS JUST GET THE
  // IMPLIMENTATION DONE
  vector<hsize_t>::iterator x_size = max_element(x_dims.begin(), x_dims.end());
  vector<hsize_t>::iterator y_size = max_element(y_dims.begin(), y_dims.end());
  // Allocate memory - resize() must be used over reserve() to ensure enough
  // buffer when we read in, _and_ that size() correctly returns the number of
  // elements in the buffer.
  f_vec->x.resize(*x_size);
  f_vec->y.resize(*y_size);
  // Now read the data into the vectors
  read_field_from_struct("f_vec", "fx_vec", f_vec->x.data());
  read_field_from_struct("f_vec", "fy_vec", f_vec->y.data());
}
