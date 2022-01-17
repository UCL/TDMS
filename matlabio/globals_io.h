//---------------------------------------------------------------------------
#ifndef globalsH
#define globalsH
//---------------------------------------------------------------------------

#include <complex>
using namespace std;

const double PI = 3.14159265358979323846;   // pi
const double eps0 = 8.8541878176e-12; // electric permittivity vacuum
const double mu0 = 4*PI*1e-7;         // magnetic permeability vacuum
const double Z0 = sqrt(mu0/eps0);     // impedance of vacuum
const complex<double> Zero(0.0,0.0);
const complex<double> Unity(1.0,0.0);
const complex<double> Two(2.0,0.0);
const complex<double> Imag(0.0,1.0);

//precission
//const double COEF_SET_ZERO = 1e-13;  //numeric_limits<double>::epsilon();
const double COEF_SET_ZERO = 0.0;  //numeric_limits<double>::epsilon();
const double PREC_D = 1e-15;  //numeric_limits<double>::epsilon();
const double ACC_D = 1e-306;  //numeric_limits<double>::denorm_min();
const double EXP_D = 706;     //maximum exponent(base-e) , 306 (base-10)
const double TOL = 1e-14;     //tolerance for calculating bessel-zeros
const double TOL2 = TOL*TOL;  //square tolerance for bessel-zeros
const int NITMX = 200;        //maximum number of itterations for bessel-zeros
const int MAXIT = 100;        //maximum number of itterations for bessel-zeros

// theta = angle of incidence between the y-axis,          xsi   = sin(theta), beta = cos(theta)
// phi   = angle between the illumination and the z-axis,  sigma = phi
// tem   = polarisation state TE=1 TM=2                    TEM   = tem
class E_in
{
  public:
  	E_in(double,double,int);
  	E_in(int,int,double,int,double);
  	E_in(int,int,double*, double*,complex<double>**,complex<double>**);
	E_in(int,double,double);
	~E_in();
  	double ret_theta();
  	double ret_phi();
	double xsi, sigma, beta, w, zr;
	int TEM, p, l, pol, Nxsi, Nsigma, type;
	complex<double> A, B;
	double *Axsi,  *Asigma;
	complex<double> **E, **H;
};


class GL
{
  public:
  	double *Xn, *Wn;
  	int N;
	GL(double x1, double x2, int n);
	GL(int n);
	~GL();
  private:
	void gauleg(double x1, double x2, int n, double *a1, double *w1);
	void gaulag(int n, double *a1, double *w1);
};

#endif
