#include "dimensions.h"
#include "tensor_init.h"
#include "matlabio.h"
#include "utils.h"

using namespace std;


void init_grid_tensors(const mxArray *ptr, SplitField &E_s, SplitField &H_s, uint8_t*** &materials){

  const char elements[][15] = {"Exy", "Exz", "Eyx", "Eyz", "Ezx", "Ezy",
                               "Hxy", "Hxz", "Hyx", "Hyz", "Hzx", "Hzy", "materials"};

  auto num_fields = mxGetNumberOfFields(ptr);
  if (num_fields != 13) {
    throw runtime_error("fdtdgrid should have 13 members, it only has " + to_string(num_fields));
  }

  for (int i = 0; i < num_fields; i++) {
    auto element = mxGetField((mxArray *) ptr, 0, elements[i]);
    string element_name = elements[i];

    if (!mxIsDouble(element) && !mxIsUint8(element)) {
      throw runtime_error("Incorrect data type in fdtdgrid. " + element_name);
    }

    auto ndims = mxGetNumberOfDimensions(element);
    if (ndims != 2 && ndims != 3){
      throw runtime_error("field matrix %s should be 2- or 3-dimensional " + element_name);
    }

    auto dims = Dimensions(element);

    if (are_equal(elements[i], "Exy")) {
      E_s.xy = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Exz")) {
      E_s.xz = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Eyx")) {
      E_s.yx = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Eyz")) {
      E_s.yz = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Ezx")) {
      E_s.zx = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Ezy")) {
      E_s.zy = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Hxy")) {
      H_s.xy = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Hxz")) {
      H_s.xz = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Hyx")) {
      H_s.yx = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Hyz")) {
      H_s.yz = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Hzx")) {
      H_s.zx = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "Hzy")) {
      H_s.zy = castMatlab3DArray(mxGetPr(element), dims[0], dims[1], dims[2]);
    } else if (are_equal(elements[i], "materials")) {

      materials = castMatlab3DArrayUint8((uint8_t *)mxGetPr(element), dims[0], dims[1], dims[2]);
      E_s.I_tot = H_s.I_tot = dims[0] - 1; // The _tot variables do NOT include the additional cell
      E_s.J_tot = H_s.J_tot = dims[1] - 1; // at the edge of the grid which is only partially used
      E_s.K_tot = H_s.K_tot = dims[2] - 1;
    } else {
      throw runtime_error("element fdtdgrid.%s not handled " + element_name);
    }
  }

  E_s.is_matlab_allocated = H_s.is_matlab_allocated = true;
}
