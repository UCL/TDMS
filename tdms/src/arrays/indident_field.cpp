#include "arrays/incident_field.h"

#include <spdlog/spdlog.h>

#include "matlabio.h"

using namespace std;

void IncidentField::set_component(Tensor3D<double> &component,
                                  const mxArray *ptr, const std::string &name) {

  if (mxIsEmpty(mxGetField(ptr, 0, name.c_str()))) {
    spdlog::info("{} not present", name);
    return;
  }

  auto element = ptr_to_nd_array_in(ptr, 3, name, "tdfield");
  auto dims = mxGetDimensions(element);
  int N = dims[0], M = dims[1], O = dims[2];
  component.initialise(cast_matlab_3D_array(mxGetPr(element), N, M, O), O, M,
                       N);

  spdlog::info("Got tdfield, dims=({},{},{})", N, M, O);
}

IncidentField::IncidentField(const mxArray *ptr) {

  assert_is_struct_with_n_fields(ptr, 2, "tdfield");
  set_component(x, ptr, "exi");
  set_component(y, ptr, "eyi");
}
