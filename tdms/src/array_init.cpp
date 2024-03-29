#include "array_init.h"

#include "dimensions.h"
#include "matlabio.h"
#include "utils.h"

using namespace std;


void init_grid_arrays(const mxArray *ptr, SplitField &E_s, SplitField &H_s,
                      uint8_t ***&materials) {

  const char elements[][15] = {"Exy", "Exz", "Eyx",      "Eyz", "Ezx",
                               "Ezy", "Hxy", "Hxz",      "Hyx", "Hyz",
                               "Hzx", "Hzy", "materials"};

  auto num_fields = mxGetNumberOfFields(ptr);
  if (num_fields != 13) {
    throw runtime_error("fdtdgrid should have 13 members, it only has " +
                        to_string(num_fields));
  }

  for (int i = 0; i < num_fields; i++) {
    auto element = mxGetField((mxArray *) ptr, 0, elements[i]);
    string element_name = elements[i];

    if (!mxIsDouble(element) && !mxIsUint8(element)) {
      throw runtime_error("Incorrect data type in fdtdgrid. " + element_name);
    }

    auto ndims = mxGetNumberOfDimensions(element);
    if (ndims != 2 && ndims != 3) {
      throw runtime_error("field matrix %s should be 2- or 3-dimensional " +
                          element_name);
    }

    auto dims = Dimensions(element);
    auto tensor =
            cast_matlab_3D_array(mxGetPr(element), dims[0], dims[1], dims[2]);

    if (are_equal(elements[i], "Exy")) {
      E_s.xy.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Exz")) {
      E_s.xz.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Eyx")) {
      E_s.yx.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Eyz")) {
      E_s.yz.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Ezx")) {
      E_s.zx.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Ezy")) {
      E_s.zy.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Hxy")) {
      H_s.xy.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Hxz")) {
      H_s.xz.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Hyx")) {
      H_s.yx.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Hyz")) {
      H_s.yz.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Hzx")) {
      H_s.zx.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "Hzy")) {
      H_s.zy.initialise_from_matlab(tensor, dims);
    } else if (are_equal(elements[i], "materials")) {
      materials = cast_matlab_3D_array((uint8_t *) mxGetPr(element), dims[0],
                                       dims[1], dims[2]);
      E_s.tot.i = H_s.tot.i =
              dims[0] -
              1;// The _tot variables do NOT include the additional cell
      E_s.tot.j = H_s.tot.j =
              dims[1] -
              1;// at the edge of the grid which is only partially used
      E_s.tot.k = H_s.tot.k = dims[2] - 1;
    } else {
      throw runtime_error("element fdtdgrid.%s not handled " + element_name);
    }
  }
}
