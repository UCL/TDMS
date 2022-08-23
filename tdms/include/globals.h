#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>

// ********
//  Macros
// ********


#define round(x)			( (int)((x)+ 0.5) ) 
#define RoundToZero(x)		( fabs(x) > 1.0e-12 ? (x) : 0.0 )


// ******************
//  Type Definitions
// ******************


typedef int*				IArray_1d;
typedef IArray_1d*			IArray_2d;
typedef IArray_2d*			IArray_3d;

typedef double*				DArray_1d;
typedef DArray_1d*			DArray_2d;
typedef DArray_2d*			DArray_3d;

typedef std::complex<double>*	CArray_1d;
typedef CArray_1d*			CArray_2d;
typedef CArray_2d*			CArray_3d;

typedef struct PlanarInterface	// Structure definition for a planar six-face interface 
{
	int I1;
	int I2;
	int J1;
	int J2;
	int K1;
	int K2;
} PlanarInterface;

typedef struct complex_vector
{
	std::complex<double> X;
	std::complex<double> Y;
	std::complex<double> Z;
} complex_vector;


// **********************
//  Enumerated constants
// **********************

enum ModeOfRun { Pass1 , Pass2 };
enum RCSType   { parallel , perpendicular };

// **************************************
//			Mathematical Constants
// **************************************

const double dcpi  = M_PI;


// **************************************
//			Physical Constants
// **************************************

const double eo	     =   8.85400e-12;				// free space electric permitivity
const double mo	     =   4.0 * dcpi * 1.0e-7;		// free space magnetic permeability
const double light_v =   1.0 / sqrt(mo*eo);			// free spave light velocity
const double Zo		 = 376.734;						// free space inpedance

const double zero			   = 1.0e-12;			// zero level for numerical comparisons
const double convergence_level = 1.0e-10;			// zero level for convergence

const auto I = std::complex<double>( 0.0 , 1.0 );
