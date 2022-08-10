#include "cmath"
#include "cstdlib"
#include "fftw3.h"
#include "globals.h"
#include "numerical_derivative.h"
#include "matlabio.h"

using namespace std;


void complex_mult_vec(fftw_complex *a, fftw_complex *b, fftw_complex *c, int len){

  for(int i=0; i<=(len-1); i++){
    // WARNING: these intermediate variables *must* be defined
    double c0 = a[i][0]*b[i][0] - a[i][1]*b[i][1];
    double c1 = a[i][0]*b[i][1] + a[i][1]*b[i][0];
    c[i][0] = c0;
    c[i][1] = c1;
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
  auto *an = (double *)malloc(sizeof(double) * N);
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
