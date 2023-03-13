#include "field.h"

#include <stdexcept>
#include <string>

#include <spdlog/spdlog.h>

#include "matlabio.h"
using std::string;

void TDFieldExporter2D::allocate(int _nI, int _nK) {

  nI = _nI;
  nK = _nK;
  mwSize dimensions[2] = {_nI, _nK};
  matlab_array = mxCreateNumericArray(2, (const mwSize *) dimensions,
                                      mxDOUBLE_CLASS, mxREAL);
  array = cast_matlab_2D_array(mxGetPr((mxArray *) matlab_array), _nI, _nK);
}

TDFieldExporter2D::~TDFieldExporter2D() { free_cast_matlab_2D_array(array); }

void TDFieldExporter2D::export_field(SplitField &F, int stride,
                                     int iteration) const {

  // need to check that enough memory was allocated before we cast!
  if ((nI < F.tot.i) || (nK < F.tot.k)) {
    throw std::runtime_error("Not enough memory to write this field! (" +
                             to_string(nI) + " , " + to_string(nK) +
                             ") but need (" + to_string(F.tot.i) + " , " +
                             to_string(F.tot.k) + ")");
  }
  // if we have enough memory, then we can write out
  int i = 0;
  while (i < F.tot.i) {
    int k = 0;
    while (k < F.tot.k) {
      array[k][i] = F.xy[k][0][i] + F.xz[k][0][i];
      k += stride;
    }
    i += stride;
  }

  // set up and write to MATLAB file
  string file_name = folder_name + "ex_" + std::to_string(iteration) + ".mat";
  spdlog::debug("Writing time-domain output to: {}", file_name);
  auto output_file = matOpen(file_name.c_str(), "w");
  matPutVariable(output_file, "ex_tdf", (mxArray *) matlab_array);
  matClose(output_file);
}
