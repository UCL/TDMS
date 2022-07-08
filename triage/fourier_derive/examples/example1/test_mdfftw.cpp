#include <fftw3.h>
#include "matrix.h"
#include "mex.h"
#include "math.h"
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include "numerical_derivative.h"
#include <complex>
using namespace std;

#include "matlabio.h"


const double dcpi  = 3.14159265358979;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{

  const mwSize *dimptr_out;
  int N;int ndims;
  double delta;
  mwSize *dims;dims = (mwSize *)malloc(3*sizeof(mwSize));
  fftw_complex *in, *out;
  fftw_plan pb, pf;
  

  
  if( nrhs != 1 ){
    fprintf(stderr,"%d\n",nrhs);
    fprintf(stderr,"Expected 1 input.");
  }
  if( nlhs != 1 ){
    fprintf(stderr,"%d\n",nlhs);
    fprintf(stderr,"Expected 1 output.");
  }
  fprintf(stderr,"position 1\n");
  double *input, *output;
  dimptr_out = mxGetDimensions(prhs[0]);
  N = mxGetNumberOfElements(prhs[0]);
  
 
  //fprintf(stderr,"N=%d\n",N);
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fprintf(stderr,"position 2\n");
  
  //pf = fftw_plan_dft_2d(dimptr_out[0], dimptr_out[1], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pf = fftw_plan_dft_2d(dimptr_out[0], dimptr_out[1], in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  fprintf(stderr,"position 3\n");
  
  ndims = 2;
  dims[0] = dimptr_out[0];
  dims[1] = dimptr_out[1]; //one for each component of field
  dims[2] = 0;
  plhs[0] = mxCreateNumericArray( ndims, (const mwSize *)dims, mxDOUBLE_CLASS, mxCOMPLEX);
  fprintf(stderr,"position 4 [%d %d %d]\n",dims[0],dims[1],N);

  for(int i=0;i<dimptr_out[0];i++)
    for(int j=0;j<dimptr_out[1];j++){
      in[j+i*dimptr_out[1]][0] = *(mxGetPr(prhs[0])+i+j*dimptr_out[0]);
      in[j+i*dimptr_out[1]][1] = *(mxGetPi(prhs[0])+i+j*dimptr_out[0]);
      
  }
  fprintf(stderr,"position 5\n");

  fftw_execute(pf);
  fprintf(stderr,"position 6\n");

  for(int i=0;i<dimptr_out[0];i++)
    for(int j=0;j<dimptr_out[1];j++){
      *(mxGetPr(plhs[0])+i+j*dimptr_out[0]) = out[j+i*dimptr_out[1]][0];
      *(mxGetPi(plhs[0])+i+j*dimptr_out[0]) = out[j+i*dimptr_out[1]][1];
      
    }
      fprintf(stderr,"position 7\n");

  //fprintf(stderr,"M02\n");
  free(dims);
  //fprintf(stderr,"M03\n");
  fftw_destroy_plan(pf);
  //fprintf(stderr,"M04\n");
  fftw_free(in);
  //fprintf(stderr,"M05\n");
  fftw_free(out);
  //fprintf(stderr,"M06\n");

}


