# include <catch2/catch_test_macros.hpp>
# include "interpolate_Efield.h"
# include "interpolate_Hfield.h"

# include <cmath>
# include <iostream>
# include <iomanip>

using namespace std;

/* Overview of Tests

We will test the performance of BLi at interpolating a known field to the centre of each Yee cell, for both the E- and H- fields.
The benchmark for success will be superior or equivalent error (Frobenius and 1D-dimension-wise norm) to that produced by MATLAB performing the same functionality.
Frobenious error is the Frobenius norm of the 3D matrix whose (i,j,k)-th entry is the error in the interpolated value at Yee cell centre (i,j,k).
The slice error, for a fixed j,k, is the norm-error of the 1D array of interpolated values along the axis (:,j,k).
The maximum of this is then the max-slice error: keeping track of this ensures us that the behaviour of BLi is consistent, and does not dramaically over-compensate in some areas and under-compensate in others.

All tests will be performed with cell sizes Dx = 0.25, Dy = 0.1, Dz = 0.05, over the range [-2,2].

NOTE: Rather than direct (<=) comparison, should impliment a relative-value check instead?
*/

inline double Frobenius(double ***M, int d1, int d2, int d3) {
    double norm_val = 0.;
    for(int i1=0; i1<d1; i1++) {
        for(int i2=0; i2<d2; i2++) {
            for(int i3=0; i3<d3; i3++) {
                norm_val += M[i1][i2][i3] * M[i1][i2][i3];
            }
        }
    }
    return sqrt(norm_val);
}
inline double Euclidean(double *v, int end, int start=0) {
    double norm_val = 0.;
    for(int i=start; i<=end; i++) {
        norm_val += v[i]*v[i];
    }
    return sqrt(norm_val);
}

inline double Ecomponent(double x, double y, double z) {
    return sin(2. * M_PI * z) * exp(-y * y) * (1. / (10 * x * x + 1.));
}
inline double Hcomponent(double x, double y, double z) {
    return cos(0.5 * M_PI * z) * exp(-y * y) * (1. / (5 * x * x + 1.));
}

// for memory allocation of 3D arrays
inline double ***allocate3dmemory(int I, int J, int K) {

    double ***p = (double ***)malloc(K * sizeof(double **));
    for (int k = 0; k < K; k++) {
        p[k] = (double **)malloc(J * sizeof(double *));
        for (int j = 0; j < J; j++) {
            p[k][j] = (double *)malloc(I * sizeof(double));
        }
    }
    return p;
}

/**
 * @brief Test the interpolation of the E-field components to the centre of the Yee cells
 *
 * Each component of the E-field will take the form
 * E{x,y,z}(xx,yy,zz) = sin(2\pi zz) * exp(-yy^2) * ( 1./ (10xx^2+1) ).
 * 
 * We test both the Fro- and slice-norm metrics, since interpolation only happens along one axis
 */
TEST_CASE("E-field interpolation check") 
{
    cout << "Beginning E-field BLi test..." << endl;
    // error tolerance, based on MATLAB performance
    double Ex_fro_tol = 2.8200485621983595e-01, Ex_ms_tol = 1.2409211493579948e-02;
    double Ey_fro_tol = 7.8295329699969822e-03, Ey_ms_tol = 7.5320765734192925e-04;
    double Ez_fro_tol = 7.5650677900775624e-03, Ez_ms_tol = 1.3131049239745484e-03;

    // additional tolerance to allow for floating-point rounding imprecisions, etc
    double acc_tol = 1e-12;

    // fake domain setup
    double cellDims[3] = {0.25, 0.1, 0.05};
    double x_lower = -2., y_lower = -2., z_lower = -2.;
    double extent_x = 4., extent_y = 4., extent_z = 4.;
    int Nx = ceil(extent_x/cellDims[0]), Ny = ceil(extent_y/cellDims[1]), Nz = ceil(extent_z/cellDims[2]);

    // setup the "split" E-field components
    double ***Exy = allocate3dmemory(Nx, Ny, Nz), ***Exz = allocate3dmemory(Nx, Ny, Nz),
           ***Eyx = allocate3dmemory(Nx, Ny, Nz), ***Eyz = allocate3dmemory(Nx, Ny, Nz),
           ***Ezx = allocate3dmemory(Nx, Ny, Nz), ***Ezy = allocate3dmemory(Nx, Ny, Nz);
    // storage for pointwise errors (-1 since we don't interpolate to cell 0's centre)
    double ***Ex_error = allocate3dmemory(Nx - 1, Ny, Nz),
           ***Ey_error = allocate3dmemory(Nx, Ny - 1, Nz),
           ***Ez_error = allocate3dmemory(Nx, Ny, Nz - 1);

    // compute the exact field and the "split field" components
    // the interpolation functions are expecting split fields, but we can bypass this by making one split field component equal to _the entire field_ value, and the other zero
    for(int ii=0; ii<Nx; ii++) {
        for(int jj=0; jj<Ny; jj++) {
            for(int kk=0; kk<Nz; kk++) {
                // these are the coordinates of Yee cell i,j,k's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5)*cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5)*cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5)*cellDims[2];

                // Initialise sample values that we will pass to the interpolation schemes:
                // in each case, set one component of the split field to be the "whole" field,
                // and the other to be 0.

                // E{x,y,z} offsets from cell centres are 0.5*D{x,y,z}
                Exy[kk][jj][ii] = Ecomponent(cell_centre[0] + 0.5*cellDims[0], cell_centre[1], cell_centre[2]);
                Exz[kk][jj][ii] = 0.;
                Eyx[kk][jj][ii] = Ecomponent(cell_centre[0], cell_centre[1] + 0.5*cellDims[1], cell_centre[2]);
                Eyz[kk][jj][ii] = 0.;
                Ezx[kk][jj][ii] = Ecomponent(cell_centre[0], cell_centre[1], cell_centre[2] + 0.5*cellDims[2]);
                Ezy[kk][jj][ii] = 0.;
            }
        }
    }

    // run interpolation functions

    // interpolate Ex
    for (int ii = 0; ii < Nx - 1; ii++)
    {
        for (int jj = 0; jj < Ny; jj++)
        {
            for (int kk = 0; kk < Nz; kk++)
            {
                // we are interpolating to the centre of cell ii+1,jj,kk
                // these are the coordinates of the Yee cell's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 1.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Ex_exact = Ecomponent(cell_centre[0], cell_centre[1], cell_centre[2]);

                // interpolate to the centre of this cell
                double Ex_interp;
                interpolateTimeDomainEx(Exy, Exz, ii + 1, jj, kk, Nx, &Ex_interp);

                // compute the errors
                Ex_error[kk][jj][ii] = Ex_interp - Ex_exact;
            }
        }
    }
    // interpolate Ey
    for (int ii = 0; ii < Nx; ii++)
    {
        for (int jj = 0; jj < Ny - 1; jj++)
        {
            for (int kk = 0; kk < Nz; kk++)
            {
                // we are interpolating to the centre of cell ii,jj+1,kk
                // these are the coordinates of the Yee cell's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 1.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Ey_exact = Ecomponent(cell_centre[0], cell_centre[1], cell_centre[2]);

                // interpolate to the centre of this cell
                double Ey_interp;
                interpolateTimeDomainEy(Eyx, Eyz, ii, jj + 1, kk, Ny, &Ey_interp);

                // compute the errors
                Ey_error[kk][jj][ii] = Ey_interp - Ey_exact;
            }
        }
    }
    // interpolate Ez
    for (int ii = 0; ii < Nx; ii++)
    {
        for (int jj = 0; jj < Ny; jj++)
        {
            for (int kk = 0; kk < Nz - 1; kk++)
            {
                // we are interpolating to the centre of cell ii,jj,kk+1
                // these are the coordinates of the Yee cell's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 1.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Ez_exact = Ecomponent(cell_centre[0], cell_centre[1], cell_centre[2]);

                // interpolate to the centre of this cell
                double Ez_interp;
                interpolateTimeDomainEz(Ezx, Ezy, ii, jj, kk + 1, Nz, &Ez_interp);

                // compute the errors
                Ez_error[kk][jj][ii] = Ez_interp - Ez_exact;
            }
        }
    }

    // can now deallocate our sample field arrays
    delete Exy, Exz, Eyx, Eyz, Ezx, Ezy;

    // compute Frobenius norms
    double Ex_fro_err = Frobenius(Ex_error, Nz, Ny, Nx - 1);
    double Ey_fro_err = Frobenius(Ey_error, Nz, Ny - 1, Nx);
    double Ez_fro_err = Frobenius(Ez_error, Nz - 1, Ny, Nx);

    // compute max-slice errors
    double Ex_ms_err = 0., Ey_ms_err = 0., Ez_ms_err = 0.;
    // Ex-slices
    for (int jj=0; jj<Ny; jj++) {
        for (int kk=0; kk<Nz; kk++) {
            // "slices" might not constitute sequential memory
            // as such, make a new array for safety
            double jk_errors[Nx-1];
            for (int ii=0; ii<Nx-1; ii++) {
                jk_errors[ii] = Ex_error[kk][jj][ii];
            }
            // compute norm-error of this slice
            double jk_slice_error = Euclidean(jk_errors, Nx-1);
            // if this exceeds the current recorded maximum error, record this
            if (jk_slice_error > Ex_ms_err) {
                Ex_ms_err = jk_slice_error;
            }
        }
    }
    // Ey-slices
    for (int ii=0; ii<Nx; ii++) {
        for (int kk=0; kk<Nz; kk++) {
            double ik_errors[Ny-1];
            for (int jj=0; jj<Ny-1; jj++) {
                ik_errors[jj] = Ey_error[kk][jj][ii];
            }
            double ik_slice_error = Euclidean(ik_errors, Ny-1);
            if (ik_slice_error > Ey_ms_err) {
                Ey_ms_err = ik_slice_error;
            }
        }
    }
    // Ez-slices
    for (int ii=0; ii<Nx; ii++) {
        for (int jj=0; jj<Ny; jj++) {
            double ij_errors[Nz-1];
            for (int kk=0; kk<Nz-1; kk++) {
                ij_errors[kk] = Ez_error[kk][jj][ii];
            }
            double ij_slice_error = Euclidean(ij_errors, Nz-1);
            if (ij_slice_error > Ez_ms_err) {
                Ez_ms_err = ij_slice_error;
            }
        }
    }

    // check Frobenius errors are acceptable
    CHECK(Ex_fro_err <= Ex_fro_tol + acc_tol);
    CHECK(Ey_fro_err <= Ey_fro_tol + acc_tol);
    CHECK(Ez_fro_err <= Ez_fro_tol + acc_tol);
    
    // check max-slice errors are acceptable
    CHECK(Ex_ms_err <= Ex_ms_tol + acc_tol);
    CHECK(Ey_ms_err <= Ey_ms_tol + acc_tol);
    CHECK(Ez_ms_err <= Ez_ms_tol + acc_tol);

    // print information to the debugger/log
    cout << " Component | Frobenius err. : (     diff     ) | Max-slice err. : (     diff     )" << endl;
    cout << fixed << scientific << setprecision(8);
    cout << "    x      | " << Ex_fro_err << " : (" << abs(Ex_fro_err - Ex_fro_tol);
    cout << ") | " << Ex_ms_err << " : (" << abs(Ex_ms_err - Ex_ms_tol) << ")" << endl;
    cout << "    y      | " << Ey_fro_err << " : (" << abs(Ey_fro_err - Ey_fro_tol);
    cout << ") | " << Ey_ms_err << " : (" << abs(Ey_ms_err - Ey_ms_tol) << ")" << endl;
    cout << "    z      | " << Ez_fro_err << " : (" << abs(Ez_fro_err - Ez_fro_tol);
    cout << ") | " << Ez_ms_err << " : (" << abs(Ez_ms_err - Ez_ms_tol) << ")" << endl;

    // memory cleanup
    delete Ex_error, Ey_error, Ez_error;
}

/**
 * @brief Test the interpolation of the H-field components to the centre of the Yee cells
 *
 * Each component of the H-field will take the form
 * H{x,y,z}(xx,yy,zz) = cos(0.5\pi zz) * exp(-yy^2) * ( 1./ (5xx^2+1) ).
 *
 * We only test Fro-norm error metrics, since interpolation must occur along two axes for each component
 */
TEST_CASE("H-field interpolation check")
{
    cout << "Beginning H-field BLi test..." <<endl;
    // error tolerance, based on MATLAB performance
    double Hx_fro_tol = 1.8564584213212786e-02;
    double Hy_fro_tol = 8.0005513762062608e-02;
    double Hz_fro_tol = 8.0012805533365733e-02;

    // additional tolerance to allow for floating-point rounding imprecisions, etc
    double acc_tol = 1e-8;

    // fake domain setup
    double cellDims[3] = {0.25, 0.1, 0.05};
    double x_lower = -2., y_lower = -2., z_lower = -2.;
    double extent_x = 4., extent_y = 4., extent_z = 4.;
    int Nx = ceil(extent_x / cellDims[0]), Ny = ceil(extent_y / cellDims[1]), Nz = ceil(extent_z / cellDims[2]);

    // setup the "split" E-field components
    double ***Hxy = allocate3dmemory(Nx, Ny, Nz), ***Hxz = allocate3dmemory(Nx, Ny, Nz),
           ***Hyx = allocate3dmemory(Nx, Ny, Nz), ***Hyz = allocate3dmemory(Nx, Ny, Nz),
           ***Hzx = allocate3dmemory(Nx, Ny, Nz), ***Hzy = allocate3dmemory(Nx, Ny, Nz);
    // storage for pointwise errors (-1 since we don't interpolate to cell 0's centre)
    double ***Hx_error = allocate3dmemory(Nx, Ny - 1, Nz - 1),
           ***Hy_error = allocate3dmemory(Nx - 1, Ny, Nz - 1),
           ***Hz_error = allocate3dmemory(Nx - 1, Ny - 1, Nz);

    // compute the exact field and the "split field" components
    // the interpolation functions are expecting split fields, but we can bypass this by making one split field component equal to _the entire field_ value, and the other zero
    for (int ii = 0; ii < Nx; ii++)
    {
        for (int jj = 0; jj < Ny; jj++)
        {
            for (int kk = 0; kk < Nz; kk++)
            {
                // these are the coordinates of Yee cell i,j,k's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // Initialise sample values that we will pass to the interpolation schemes:
                // in each case, set one component of the split field to be the "whole" field,
                // and the other to be 0.

                // H{x,y,z} offsets from cell centres are 0.5 * D{!{x,y,z}} IE, 0.5 away from the centre in the two directions that are _not_ the field component
                Hxy[kk][jj][ii] = Hcomponent(cell_centre[0],
                                             cell_centre[1] + 0.5 * cellDims[1],
                                             cell_centre[2] + 0.5 * cellDims[2]);
                Hxz[kk][jj][ii] = 0.;
                Hyx[kk][jj][ii] = Hcomponent(cell_centre[0] + 0.5 * cellDims[0],
                                             cell_centre[1],
                                             cell_centre[2] + 0.5 * cellDims[2]);
                Hyz[kk][jj][ii] = 0.;
                Hzx[kk][jj][ii] = Hcomponent(cell_centre[0] + 0.5 * cellDims[0],
                                             cell_centre[1] + 0.5 * cellDims[1],
                                             cell_centre[2]);
                Hzy[kk][jj][ii] = 0.;
            }
        }
    }

    // run interpolation functions
    
    // interpolate Hx
    for (int ii = 0; ii < Nx; ii++)
    {
        for (int jj = 0; jj < Ny - 1; jj++)
        {
            for (int kk = 0; kk < Nz - 1; kk++)
            {
                // we are interpolating to the centre of cell ii,jj+1,kk+1
                // these are the coordinates of Yee cell centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 1.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 1.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Hx_exact = Hcomponent(cell_centre[0], cell_centre[1], cell_centre[2]);

                // interpolate to the centre of this cell
                double Hx_interp;
                interpolateTimeDomainHx(Hxy, Hxz, ii, jj+1, kk+1, Ny, Nz, &Hx_interp);

                // compute the errors
                Hx_error[kk][jj][ii] = Hx_interp - Hx_exact;
            }
        }
    }
    // interpolate Hy
    for (int ii = 0; ii < Nx - 1; ii++)
    {
        for (int jj = 0; jj < Ny; jj++)
        {
            for (int kk = 0; kk < Nz - 1; kk++)
            {
                // we are interpolating to the centre of cell ii,jj+1,kk+1
                // these are the coordinates of Yee cell centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 1.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 1.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Hy_exact = Hcomponent(cell_centre[0], cell_centre[1], cell_centre[2]);

                // interpolate to the centre of this cell
                double Hy_interp;
                interpolateTimeDomainHy(Hyx, Hyz, ii + 1, jj, kk + 1, Nx, Nz, &Hy_interp);

                // compute the errors
                Hy_error[kk][jj][ii] = Hy_interp - Hy_exact;
            }
        }
    }
    // interpolate Hz
    for (int ii = 0; ii < Nx - 1; ii++)
    {
        for (int jj = 0; jj < Ny - 1; jj++)
        {
            for (int kk = 0; kk < Nz; kk++)
            {
                // we are interpolating to the centre of cell ii+1,jj+1,kk
                // these are the coordinates of Yee cell centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 1.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 1.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Hz_exact = Hcomponent(cell_centre[0], cell_centre[1], cell_centre[2]);

                // interpolate to the centre of this cell
                double Hz_interp;
                interpolateTimeDomainHz(Hzx, Hzy, ii + 1, jj + 1, kk, Nx, Ny, &Hz_interp);

                // compute the errors
                Hz_error[kk][jj][ii] = Hz_interp - Hz_exact;
            }
        }
    }

    // can now deallocate our sample field arrays
    delete Hxy, Hxz, Hyx, Hyz, Hzx, Hzy;
    
    // compute Frobenius norms
    double Hx_fro_err = Frobenius(Hx_error, Nz - 1, Ny - 1, Nx);
    double Hy_fro_err = Frobenius(Hy_error, Nz - 1, Ny, Nx - 1);
    double Hz_fro_err = Frobenius(Hz_error, Nz, Ny - 1, Nx - 1);

    // check Frobenius errors are acceptable
    CHECK(Hx_fro_err <= Hx_fro_tol + acc_tol);
    CHECK(Hy_fro_err <= Hy_fro_tol + acc_tol);
    CHECK(Hz_fro_err <= Hz_fro_tol + acc_tol);

    // print information to the debugger/log
    cout << " Component | Frobenius err. : (     diff     )" << endl;
    cout << fixed << scientific << setprecision(8);
    cout << "    x      | " << Hx_fro_err << " : (" << abs(Hx_fro_err - Hx_fro_tol) << ")" << endl;
    cout << "    y      | " << Hy_fro_err << " : (" << abs(Hy_fro_err - Hy_fro_tol) << ")" << endl;
    cout << "    z      | " << Hz_fro_err << " : (" << abs(Hz_fro_err - Hz_fro_tol) << ")" << endl;

    // memory cleanup
    delete Hx_error, Hy_error, Hz_error;
}
