#include "matlabio.h"
#include "field.h"


void TDFieldExporter2D::allocate(int I, int J) {

  mwSize dims[2] = {I, J};
  matlab_array = mxCreateNumericArray(2, (const mwSize *)dims, mxDOUBLE_CLASS, mxREAL);
  array = castMatlab2DArray(mxGetPr((mxArray *) matlab_array), I, J);
}

TDFieldExporter2D::~TDFieldExporter2D() {
  freeCastMatlab2DArray(array);
}


void TDFieldExporter2D::export_field(SplitField& F, int stride, int iteration) const{

  int ic = 0;

  for (int i = 0; i < F.I_tot; i++) {
    int kc = 0;
    if ((i % stride) == 0) {
      for (int k = 0; k < F.K_tot; k++)
        if ((k % stride) == 0) {
          array[kc++][ic] = F.xy[k][0][i] + F.xz[k][0][i];
        }
      ic++;
    }
  }

  char toutputfilename[512];
  sprintf(toutputfilename, "%s/ex_%06d.mat", folder_name, iteration);
  fprintf(stderr, "time domain output: %s\n", toutputfilename);

  auto toutfile = matOpen(toutputfilename, "w");
  matPutVariable(toutfile, "ex_tdf", (mxArray *) matlab_array);
  matClose(toutfile);
}
