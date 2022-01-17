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

void complex_mult_vec( fftw_complex *a, fftw_complex *b, fftw_complex *c, int len){
  double c0, c1;
  for(int i=0;i<=(len-1);i++){
    c0 = a[i][0]*b[i][0] - a[i][1]*b[i][1];
    c1 = a[i][0]*b[i][1] + a[i][1]*b[i][0];
    c[i][0] = c0;
    c[i][1] = c1;
    //c[i][0] = a[i][0]*b[i][0] - a[i][1]*b[i][1];
    //c[i][1] = a[i][0]*b[i][1] + a[i][1]*b[i][0];
     
  }

}

/*
  Initialise the coefficients required to simultaneously perform 
  differentiation and shifting, by amount delta, using a forward 
  and backward FFT.
  
  delta is a fraction of the spatial step

  Dk is the coefficients

  N is the number of elements in Dk

*/
void init_diff_shift_op( double delta, fftw_complex *Dk, int N){
  //fprintf(stdout,"init_diff_shift_op, delta=%e\n",delta);
  //define an
  double *an = (double *)malloc(sizeof(double) * N);
  if( (N % 2) == 0){//even case
    for(int i=0;i<=(N/2);i++)
      an[i] = (double) i;
    //    an[N/2] = 0.;
    for(int i=N/2+1;i<=(N-1);i++)
      an[i] = (double) (i-N);
  }
  else{//odd case
    for(int i=0;i<=(N-1)/2;i++)
      an[i] = (double) i;
    for(int i=(N+1)/2;i<=(N-1);i++)
      an[i] = (double) (i-N);
  }

  for(int i=0;i<=(N-1);i++){
    Dk[i][0] = -1.*sin(an[i]*2.*dcpi*delta/( (double)N ))*an[i]*2.*dcpi;
    Dk[i][1] =  cos(an[i]*2.*dcpi*delta/( (double)N ))*an[i]*2.*dcpi;
  }
  if( (N % 2) == 0){
    Dk[N/2][0] = -1.*sin(an[N/2]*2.*dcpi*delta/( (double)N ))*an[N/2]*2.*dcpi;
    Dk[N/2][1] = 0.;
  }
  /*
  fprintf(stdout,"delta: %e\n",delta);
  for(int i=0;i<=(N-1);i++){
    fprintf(stdout,"%e %e %e\n",an[i],Dk[i][0],Dk[i][1]);
  }
  */
  free(an);
}

/* 
   Calculate the first derivative of a sampled function

   in_pb_pf: the buffer containing the data to be differentiated. 
   Note that it must be the buffer which is the input for both plans 
   pf and pb

   out_pb_pf: the buffer which will contain the computed derivative. 
   It must be the output for both plans pf and pb

   pf, pb: plans for forward and backward FFT, respectively.

   N: number of elements in buffers
 */

void first_derivative( fftw_complex *in_pb_pf, fftw_complex *out_pb_pf,
		       fftw_complex *Dk, int N, fftw_plan pf, fftw_plan pb){
  
  //fprintf(stderr,"01\n");
  fftw_execute(pf);
  //fprintf(stderr,"02\n");
  complex_mult_vec(out_pb_pf, Dk, in_pb_pf, N);
  //fprintf(stderr,"03\n");
  fftw_execute(pb);
  //fprintf(stderr,"04\n");
  for(int i=0;i<=(N-1);i++){
    out_pb_pf[i][0] = out_pb_pf[i][0]/((double) N);out_pb_pf[i][1] = out_pb_pf[i][1]/((double) N);
  }
  //fprintf(stderr,"05\n");

}

/*
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{

  const mwSize *dimptr_out;
  int N;int ndims;
  double delta;
  mwSize *dims;dims = (mwSize *)malloc(3*sizeof(mwSize));
  fftw_complex *in, *out, *zk;
  fftw_plan pb, pf;
  

  
  if( nrhs != 2 ){
    fprintf(stderr,"%d\n",nrhs);
    fprintf(stderr,"Expected 1 input.");
  }
  if( nlhs != 1 ){
    fprintf(stderr,"%d\n",nlhs);
    fprintf(stderr,"Expected 1 output.");
  }

  double *input, *output;
  dimptr_out = mxGetDimensions(prhs[0]);
  if(dimptr_out[0]==1)
    N = dimptr_out[1];
  else
    N = dimptr_out[0];
  fprintf(stderr,"N=%d\n",N);
  delta = *mxGetPr(prhs[1]);
  fprintf(stderr,"delta=%e\n",delta);
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  zk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  pf = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pb = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  ndims = 2;
  dims[0] = 1;
  dims[1] = N; //one for each component of field
  dims[2] = 0;
  plhs[0] = mxCreateNumericArray( ndims, (const mwSize *)dims, mxDOUBLE_CLASS, mxREAL);
  
  for(int i=0;i<N;i++){
    in[i][0] = *(i+mxGetPr(prhs[0])); 
    in[i][1] = 0.;
  }
  init_diff_shift_op( delta, zk, N);
  first_derivative( in, out, zk, N, pf, pb);
  //fprintf(stderr,"M01\n");
  for(int i=0;i<N;i++){
    *(i+mxGetPr(plhs[0])) = out[i][0];
  }
  //fprintf(stderr,"M02\n");
  free(dims);
  //fprintf(stderr,"M03\n");
  fftw_destroy_plan(pf);fftw_destroy_plan(pb);
  //fprintf(stderr,"M04\n");
  fftw_free(in);
  //fprintf(stderr,"M05\n");
  fftw_free(out);
  //fprintf(stderr,"M06\n");
  fftw_free(zk);
  //fprintf(stderr,"M07\n");
}


*/
