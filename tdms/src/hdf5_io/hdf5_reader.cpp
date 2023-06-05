#include "hdf5_io/hdf5_reader.h"
#include "hdf5_io/hdf5_dimension.h"

#include <algorithm>
#include <stdexcept>

#include <spdlog/spdlog.h>

using namespace std;

void HDF5Reader::read(const string &plane, InterfaceComponent &ic) const {
  string path_to_data = "interface/" + plane;
  if (path_exists("interface/" + plane, H5I_DATASET, true)) {
    // Read the InterfaceComponent in as a 2-element double array
    double read_buffer[2];
    read("interface/" + plane, read_buffer);
    // The index that is read in should have 1 subtracted from it, to account
    // for MATLAB indexing
    ic.index = max((int) read_buffer[0] - 1, 0);
    // The apply flag should be cast from the double that is read in
    ic.apply = (bool) read_buffer[1];
  }
}

void HDF5Reader::read(FrequencyVectors &f_vec) const {
  // Allocate memory in f_vec
  H5Dimension x_dims = shape_of("f_vec/fx_vec");
  H5Dimension y_dims = shape_of("f_vec/fy_vec");
  // Check that we have one dimensional arrays
  if (!x_dims.is_1D() || !y_dims.is_1D()) {
    throw runtime_error("f_vec members are not 1D arrays!");
  }
  // Now read the data into the vectors
  read("f_vec/fx_vec", f_vec.x);
  read("f_vec/fy_vec", f_vec.y);
}

void HDF5Reader::read(Cuboid &cube) const {
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
  // USE THE utils::vector functions here when possible!
  read(cuboid_dataset, intermediate_buffer);
  for (int i = 0; i < 6; i++) {
    cube.array[i] = (int) intermediate_buffer[i] - 1;
  }
}

void HDF5Reader::read(DispersiveMultiLayer &dml) const {
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
  read(group_name + "/alpha", dml.alpha);
  read(group_name + "/beta", dml.beta);
  read(group_name + "/gamma", dml.gamma);
  read(group_name + "/kappa_x", dml.kappa.x);
  read(group_name + "/kappa_y", dml.kappa.y);
  read(group_name + "/kappa_z", dml.kappa.z);
  read(group_name + "/sigma_x", dml.sigma.x);
  read(group_name + "/sigma_y", dml.sigma.y);
  read(group_name + "/sigma_z", dml.sigma.z);
}

void HDF5Reader::read(const string &group_name, const string &name_prefix,
                      XYZVector &v) const {
  for (const char &component : {'x', 'y', 'z'}) {
    read(group_name + "/" + name_prefix + component, v[component]);
  }
}

void HDF5Reader::read(CMaterial &c_material, const string &group_name) const {
  // We should expect a group with 6 members
  H5::Group group = file_->openGroup(group_name);
  int n_members = group.getNumObjs();
  if (n_members != 9) {
    throw runtime_error("CMaterial should have 9 members, but " +
                        to_string(n_members) + " were found");
  }
  group.close();

  read(group_name, "Ca", c_material.a);
  read(group_name, "Cb", c_material.b);
  read(group_name, "Cc", c_material.c);
}

void HDF5Reader::read(CCollection &c_collection,
                      const std::string &group_name) const {
  // We should expect a group with either 6 or 9 members
  H5::Group group = file_->openGroup(group_name);
  int n_members = group.getNumObjs();

  if (n_members != 6 && n_members != 9) {
    throw runtime_error("CCollection should have 6 or 9 members, but " +
                        to_string(n_members) + " were found");
  }
  group.close();

  // Initialise the a and b arrays
  read(group_name, "Ca", c_collection.a);
  read(group_name, "Cb", c_collection.b);

  // Initialise the c arrays, if provided
  // Additionally flag the material as a dispersive multilayer
  if (n_members == 9) {
    c_collection.is_disp_ml = true;
    read(group_name, "Cc", c_collection.c);
  }
}
