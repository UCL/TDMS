#include "numerical_derivative.h"

#include <cmath>
#include <cstdlib>

#include <fftw3.h>
#include <spdlog/spdlog.h>

#include "globals.h"
#include "matlabio.h"

using namespace std;
using namespace tdms_math_constants;

// fftw_complex is typdef to a double[2] - first element is Re, second Im.
const int REAL=0, IMAG=1;

void complex_mult_vec(fftw_complex *a, fftw_complex *b, fftw_complex *c, int len){

  for(int i=0; i<=(len-1); i++){
    // WARNING: these intermediate variables *must* be defined
    double c0 = a[i][REAL]*b[i][REAL] - a[i][IMAG]*b[i][IMAG];
    double c1 = a[i][REAL]*b[i][IMAG] + a[i][IMAG]*b[i][REAL];
    c[i][REAL] = c0;
    c[i][IMAG] = c1;
  }
}

void init_diff_shift_op(double delta, fftw_complex *Dk, int N){

  spdlog::debug("init_diff_shift_op, delta={}", delta);

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
    Dk[i][REAL] = -1.*sin(an[i]*2.*DCPI*delta/( (double)N ))*an[i]*2.*DCPI;
    Dk[i][IMAG] =  cos(an[i]*2.*DCPI*delta/( (double)N ))*an[i]*2.*DCPI;
  }
  if( (N % 2) == 0){
    Dk[N/2][REAL] = -1.*sin(an[N/2]*2.*DCPI*delta/( (double)N ))*an[N/2]*2.*DCPI;
    Dk[N/2][IMAG] = 0.;
  }
  spdlog::trace("delta: {}\n", delta);
  for(int i=0;i<=(N-1);i++){
    spdlog::trace("{} {} {}",an[i],Dk[i][0],Dk[i][1]);
  }
  free(an);
}

void first_derivative( fftw_complex *in_pb_pf, fftw_complex *out_pb_pf,
		       fftw_complex *Dk, int N, fftw_plan pf, fftw_plan pb){
  /*
   *  1. Fourier transform the data in in_pb_pf, placing the result into out_pb_pf
   *  2. Multiply the result of the FT (out_pb_pf) by the Dk coefficients, placing
   *     the result into in_pb_pf.
   *  3. Fourier transform back.
   *
   *  Due to how we've organised things, we can use the same fftw_scheme as in the
   *  first step, providing the data in_pb_pf and getting the output in out_pb_pf
   */
  fftw_execute(pf);
  complex_mult_vec(out_pb_pf, Dk, in_pb_pf, N);
  fftw_execute(pb);
  for(int i=0;i<=(N-1);i++){
    out_pb_pf[i][REAL] = out_pb_pf[i][REAL]/((double) N);
    out_pb_pf[i][IMAG] = out_pb_pf[i][IMAG]/((double) N);
  }
}
