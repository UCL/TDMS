#include "matlabio.h"
#include "field.h"


void TDFieldExporter2D::allocate(int nI, int nJ) {

  mwSize dims[2] = {nI, nJ};
  matlab_array = mxCreateNumericArray(2, (const mwSize *)dims, mxDOUBLE_CLASS, mxREAL);
  array = castMatlab2DArray(mxGetPr((mxArray *) matlab_array), nI, nJ);
}

TDFieldExporter2D::~TDFieldExporter2D() {
  freeCastMatlab2DArray(array);
}


void TDFieldExporter2D::export_field(SplitField& F, int stride, int iteration) const{

  int i = 0;
  while (i < F.I_tot) {
      int k = 0;
      while (k < F.K_tot){
        array[k][i] = F.xy[k][0][i] + F.xz[k][0][i];
        k += stride;
      }
      i += stride;
    }

  char toutputfilename[512];
  sprintf(toutputfilename, "%s/ex_%06d.mat", folder_name, iteration);
  fprintf(stderr, "time domain output: %s\n", toutputfilename);

  auto toutfile = matOpen(toutputfilename, "w");
  matPutVariable(toutfile, "ex_tdf", (mxArray *) matlab_array);
  matClose(toutfile);
}
