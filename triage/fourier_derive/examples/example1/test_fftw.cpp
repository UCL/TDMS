#include <stdio.h>
#include <math.h>
#include <string.h> 
#include <stdlib.h> 
//#include <complex.h>
#include <fftw3.h>
#include <complex>
//using namespace std;
void complex_mult_vec( fftw_complex *a, fftw_complex *b, fftw_complex *c, int len);
void first_derivative( fftw_complex *in_pb_pf, fftw_complex *out_pb_pf,
		       fftw_complex *Dk, int N, fftw_plan pf, fftw_plan pb);
void init_diff_shift_op( double delta, fftw_complex *Dk, int N);
//const complex<double> I = complex<double>( 0.0 , 1.0 );
const double dcpi  = 3.14159265358979;

int main(int nargs,char *argv[], char **envp){
  fftw_complex *in, *out, *yk, *xk, *zk, *qk;
  fftw_plan pb, pf;
  int N = 32;
  double delta = 0.5;
  //complex<double> I = complex<double>( 0.0 , 1.0 );
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*N));
  yk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*N));
  xk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*N));
  zk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*N));
  qk = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*N));
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*N));
  pf = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  pb = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  
  for(int i=0;i<=(N-1);i++){
    in[i][0] = (cos(2.*dcpi*( (double) i)/((double) N)) + sin(2.*dcpi*2.*( (double) i)/((double) N)));//*exp(-1.*pow( 10.*((double) (i-16.))/((double) N),2.));
    in[i][1] = 0.;
  }
  for(int i=0;i<=(N-1);i++){
    fprintf(stdout,"%e %e %e\n",( (double) i)/((double) N),in[i][0],in[i][1]); 
  }
  //  fprintf(stdout,"even?: %d\n", (N % 2) == 0);
  /*
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
    yk[i][0] = 0.;
    yk[i][1] = 2.*dcpi*an[i];
  }
  for(int i=0;i<=(N-1);i++){
    xk[i][0] = cos(an[i]*2.*dcpi*delta);
    xk[i][1] = sin(an[i]*2.*dcpi*delta);
  }
  for(int i=0;i<=(N-1);i++){
    zk[i][0] = yk[i][0]*xk[i][0] - yk[i][1]*xk[i][1];
    zk[i][1] = yk[i][0]*xk[i][1] + yk[i][1]*xk[i][0];
  }
  for(int i=0;i<=(N-1);i++){
    qk[i][0] = -sin(an[i]*2.*dcpi*delta/( (double)N ))*an[i]*2.*dcpi;
    qk[i][1] =  cos(an[i]*2.*dcpi*delta/( (double)N ))*an[i]*2.*dcpi;
  }
  if( (N % 2) == 0){
    qk[N/2][0] = 1.*sin(an[N/2]*2.*dcpi*delta/( (double)N ))*an[N/2]*2.*dcpi;
    qk[N/2][1] = 0.;
    fprintf(stderr,"qk[N/2] = %e\n",qk[N/2][0]);
  }
  

  fftw_execute(pf); //repeat as needed 
  fprintf(stderr,"G[N/2] = %e+i%e\n",out[N/2][0],out[N/2][1]);
  /*
  /*
    for(int i=0;i<=(N-1);i++){
    in[i][0] = qk[i][0]*out[i][0] - qk[i][1]*out[i][1];
    in[i][1] = qk[i][0]*out[i][1] + qk[i][1]*out[i][0];
  }
  */
  /*
  complex_mult_vec(out,qk,in,N);
  fftw_execute(pb);
  for(int i=0;i<=(N-1);i++){
    out[i][0] = out[i][0]/((double) N);out[i][1] = out[i][1]/((double) N);
  } 
  */
  /*
  for(int i=0;N-1;i++){
    in[i] = 0.;
  */
  init_diff_shift_op( delta, zk, N);
  first_derivative( in, out, zk, N, pf, pb);

    
  for(int i=0;i<=(N-1);i++){
    fprintf(stdout,"%e %e %e\n",( (double) i + delta)/((double) N),out[i][0],out[i][1]); 
  }
  for(int i=0;i<=(N-1);i++){
    fprintf(stdout,"%e %e %e\n",( (double) i + delta)/((double) N),zk[i][0],zk[i][1]); 
  }
  fftw_destroy_plan(pf);fftw_destroy_plan(pb);
  fftw_free(in); fftw_free(out);fftw_free(yk);fftw_free(xk);fftw_free(zk);fftw_free(qk);
  //free(an);

}

//c = b * a
void complex_mult_vec( fftw_complex *a, fftw_complex *b, fftw_complex *c, int len){
  for(int i=0;i<=(len-1);i++){
    c[i][0] = a[i][0]*b[i][0] - a[i][1]*b[i][1];
    c[i][1] = a[i][0]*b[i][1] + a[i][1]*b[i][0];
     
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
    Dk[i][0] = -sin(an[i]*2.*dcpi*delta/( (double)N ))*an[i]*2.*dcpi;
    Dk[i][1] =  cos(an[i]*2.*dcpi*delta/( (double)N ))*an[i]*2.*dcpi;
  }
  if( (N % 2) == 0){
    Dk[N/2][0] = -1.*sin(an[N/2]*2.*dcpi*delta/( (double)N ))*an[N/2]*2.*dcpi;
    Dk[N/2][1] = 0.;
  }

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
  
  fftw_execute(pf);
  complex_mult_vec(out_pb_pf, Dk, in_pb_pf, N);
  fftw_execute(pb);
  for(int i=0;i<=(N-1);i++){
    out_pb_pf[i][0] = out_pb_pf[i][0]/((double) N);out_pb_pf[i][1] = out_pb_pf[i][1]/((double) N);
  }
  

}
