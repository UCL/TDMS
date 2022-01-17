#include "math.h"
#include "mex.h"
#include <complex>
#include "matrix.h"

using namespace std;
#include "mesh_base.h"
#include "matlabio.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[]){
  int I0, I1, J0, J1, K,K0,K1, counter = 0,coordmap[3];//triangulate plane
  int *p1, *p2, *p3, v1[3], v2[3], *cross, dims[2];//testing cross product

  
  if(nrhs==6){//triangulate plane
    
    if(nrhs != 6)
      mexErrMsgTxt("Incorrect number of input parameters");
    
    
    I0 = ((int) *mxGetPr(prhs[counter++]));
    I1 = ((int) *mxGetPr(prhs[counter++]));
    J0 = ((int) *mxGetPr(prhs[counter++]));
    J1 = ((int) *mxGetPr(prhs[counter++]));
    K0 = ((int) *mxGetPr(prhs[counter++]));
    K1 = ((int) *mxGetPr(prhs[counter++]));
    
    if( nlhs != 2)
      mexErrMsgTxt("Must have 2 output argument");
    conciseTriangulateCuboid(I0, I1, J0, J1, K0, K1, (mxArray **)&plhs[0], (mxArray **)&plhs[1]);
  }
  else if(nrhs==3){//testing cross product
    
   if( nlhs != 1)
      mexErrMsgTxt("Must have 1 output argument");
   
   p1 = (int *)mxGetPr(prhs[0]);
   p2 = (int *)mxGetPr(prhs[1]);
   p3 = (int *)mxGetPr(prhs[2]);

   dims[0] = 1;
   dims[1] = 3;

   plhs[0] = mxCreateNumericArray( 2, dims, mxINT32_CLASS, mxREAL);
   cross = (int *)mxGetPr(plhs[0]);
   
   pointsToVector(p1, p2, v1);
   pointsToVector(p1, p3, v2);

   crossProduct( v1, v2, cross);
   
  }
  else{
    mexErrMsgTxt("Incorrect number of input parameters");
  }


}
