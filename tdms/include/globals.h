/**
 * @file globals.h
 * @brief Type definitions and global constants.
 */
#pragma once

#include <complex>

// ******************
//  Type Definitions
// ******************

typedef int *IArray_1d;
typedef IArray_1d *IArray_2d;
typedef IArray_2d *IArray_3d;

typedef double *DArray_1d;
typedef DArray_1d *DArray_2d;
typedef DArray_2d *DArray_3d;

typedef std::complex<double> *CArray_1d;
typedef CArray_1d *CArray_2d;
typedef CArray_2d *CArray_3d;

typedef struct PlanarInterface// Structure definition for a planar six-face
                              // interface
{
  int I1;
  int I2;
  int J1;
  int J2;
  int K1;
  int K2;
} PlanarInterface;

typedef struct complex_vector {
  std::complex<double> X;
  std::complex<double> Y;
  std::complex<double> Z;
} complex_vector;

enum AxialDirection { X = 'x', Y = 'y', Z = 'z' };

/**
 * Enum defining a mapping to integers used in the MATLAB initialisation
 */
enum FieldComponents { Ex = 1, Ey, Ez, Hx, Hy, Hz };

// **********************
//  Enumerated constants
// **********************

enum ModeOfRun { Pass1, Pass2 };
enum RCSType { parallel, perpendicular };
enum SolverMethod { PseudoSpectral, FiniteDifference };

//! The interpolation methods that can be used to extract field values at Yee
//! cell centres
enum PreferredInterpolationMethods { BandLimited, Cubic };

// **************************************
//			Mathematical Constants
// **************************************

namespace tdms_math_constants {
const double DCPI = 3.14159265358979323846;// Pi
const std::complex<double> IMAGINARY_UNIT =
        std::complex<double>(0.0, 1.0);// Imaginary unit
}// namespace tdms_math_constants

// **************************************
//			Physical Constants
// **************************************

namespace tdms_phys_constants {
const double EPSILON0 = 8.85400e-12;// free space electric permitivity
const double MU0 = 4.0 * tdms_math_constants::DCPI *
                   1.0e-7;// free space magnetic permeability
const double LIGHT_V = 1.0 / sqrt(EPSILON0 * MU0);// free space light velocity
const double Z0 = 376.734;                        // free space inpedance
}// namespace tdms_phys_constants
