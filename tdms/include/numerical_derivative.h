#include <fftw3.h>
void complex_mult_vec( fftw_complex *a, fftw_complex *b, fftw_complex *c, int len);
void init_diff_shift_op( double delta, fftw_complex *Dk, int N);
void first_derivative( fftw_complex *in_pb_pf, fftw_complex *out_pb_pf,
		       fftw_complex *Dk, int N, fftw_plan pf, fftw_plan pb);
