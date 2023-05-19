#include "hdf5_io/hdf5_dimension.h"

#include <algorithm>

using namespace std;

H5Dimension::H5Dimension(const H5::DataSpace &data_space) {
  // Fetch rank to declare the vector size
  int rank = data_space.getSimpleExtentNdims();
  // Resize to the correct rank then populate entries with the data
  resize(rank);
  data_space.getSimpleExtentDims(data(), nullptr);
}

bool H5Dimension::is_1D() const {
  int n_non_trivial_dimensions = 0;
  if (size() != 1) {
    for (hsize_t dim : *this) {
      if (dim > 1) { n_non_trivial_dimensions++; }
    }
  }
  return n_non_trivial_dimensions <= 1;
}
