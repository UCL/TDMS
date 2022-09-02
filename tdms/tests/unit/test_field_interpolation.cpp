# include <catch2/catch_test_macros.hpp>
# include "interpolate_Efield.h"
# include "interpolate_Hfield.h"

# include <cmath>
# include <iostream>

using namespace std;

/* OVERVIEW OF TESTS

Assume there is no illumination origin offset.
Let Dx, Dy, Dz be the dimensions of the Yee cell, so cell (i,j,k) has centre at position O_{ijk} = (i Dx, j Dy, k Dz).
The components associated to the cell (i,j,k) are located at O_{ijk} + vOffset, where vOffset is given by

Component           | vOffset ./ (Dx, Dy, Dz)
    Ex              | (0.5, 0.0, 0.0)
    Ey              | (0.0, 0.5, 0.0)
    Ez              | (0.0, 0.0, 0.5)
    Hx              | (0.5, 0.0, 0.5)
    Hy              | (0.5, 0.0, 0.5)
    Hz              | (0.5, 0.5, 0.0)

Let Nx, Ny, Nz be integers. Set:
Dx = 1/(Nx-1), Dy = 1/(Ny-1), Dz = 1/(Nz-1),

and let x,y,z be arrays of 2N_{x,y,x}+1 elements such that:
x[i] = i*Dx/2, y[j] = j*Dy/2, z[k] = k*Dz/2.

Using the slice notation start:step:stop (inclusive), this means that:
{x,y,z}[0:2:end] triples give the coordinates of the Yee cell centres
{x,y,z}[1:2:end] give the coordinates of the samples associated to the Yee cells.

Now use the notation F[i,j,k] =  F(x[i],y[j],z[k]) for some known vector field F.
For ii,jj,kk from 0 to N_{x,y,z} inclusive, we have that:
{E,H}{x,y,z}[{2*ii, 2*jj, 2*kk] is the value of {E,H}{x,y,z} at the centre of Yee cell (i,j,k)
EX_{ii,jj,kk} := Ex[2*ii+1, 2*jj, 2*kk] is the value of Ex associated to Yee cell (i,j,k)
EY_{ii,jj,kk} := Ey[2*ii, 2*jj+1, 2*kk] is the value of Ey associated to Yee cell (i,j,k)
EZ_{ii,jj,kk} := Ez[2*ii, 2*jj, 2*kk+1] is the value of Ez associated to Yee cell (i,j,k)
HX_{ii,jj,kk} := Hx[2*ii+2, 2*jj, 2*kk+1] is the value of Hx associated to Yee cell (i,j,k)
HY_{ii,jj,kk} := Hy[2*ii+1, 2*jj, 2*kk+1] is the value of Hy associated to Yee cell (i,j,k)
HZ_{ii,jj,kk} := Hz[2*ii+1, 2*jj+1, 2*kk] is the value of Hx associated to Yee cell (i,j,k)

As such, we can perform interpolation on fields that we define, and compare to the exact values of the field at those given points. */

/*
We will define the E-field to be
 * Ex(x,y,z) = cos(2\pi x)sin(2\pi y)sin(2\pi z),
 * Ey(x,y,z) = sin(2\pi x)cos(2\pi y)sin(2\pi z),
 * Ez(x,y,z) = sin(2\pi x)sin(2\pi y)cos(2\pi z),
and the H-field by
 * Hx(x,y,z) = sin(2\pi x)cos(2\pi y)cos(2\pi z),
 * Hy(x,y,z) = cos(2\pi x)sin(2\pi y)cos(2\pi z),
 * Hz(x,y,z) = cos(2\pi x)cos(2\pi y)sin(2\pi z),
*/
inline double s2pi(double x) {
    return sin(2.*M_PI*x);
}
inline double c2pi(double x) {
    return cos(2.*M_PI*x);
}
inline double Ex(double x, double y, double z) {
    return c2pi(x) * s2pi(y) * s2pi(z);
}
inline double Ey(double x, double y, double z) {
    return s2pi(x) * c2pi(y) * s2pi(z);
}
inline double Ez(double x, double y, double z) {
    return s2pi(x) * s2pi(y) * c2pi(z);
}
inline double Hx(double x, double y, double z) {
    return s2pi(x) * c2pi(y) * c2pi(z);
}
inline double Hy(double x, double y, double z) {
    return c2pi(x) * s2pi(y) * c2pi(z);
}
inline double Hz(double x, double y, double z) {
    return c2pi(x) * c2pi(y) * s2pi(z);
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
 */
TEST_CASE("E-field interpolation check") {

    // error tolerance
    // this needs to be set based on some kind of reference value.
    // Implementing this scheme in MATLAB (for the fields above) yields an error of
    // Ex : 6.30866994e-02 | Ey : 1.26746769e-01 | Ez : 2.55475206e-01
    double Ex_tol = 6.30866994e-02, Ey_tol = 1.26746769e-01, Ez_tol = 2.55475206e-01;

    // fake domain setup
    int Nx = 100, Ny = 50, Nz = 25;
    double Dx = 1./(double)Nx, Dy = 1./(double)Ny, Dz = 1./(double)Nz;

    // setup the gridpoints
    double x[2*Nx+1], y[2*Ny+1], z[2*Nz+1];
    for(int i=0; i<2*Nx+1; i++) {x[i] = i*Dx / 2.;}
    for(int j=0; j<2*Ny+1; j++) {y[j] = j*Dy / 2.;}
    for(int k=0; k<2*Nz+1; k++) {z[k] = k*Dz / 2.;}

    // setup the "split" E-field components
    double ***Exy = allocate3dmemory(Nx, Ny, Nz), ***Exz = allocate3dmemory(Nx, Ny, Nz),
           ***Eyx = allocate3dmemory(Nx, Ny, Nz), ***Eyz = allocate3dmemory(Nx, Ny, Nz),
           ***Ezx = allocate3dmemory(Nx, Ny, Nz), ***Ezy = allocate3dmemory(Nx, Ny, Nz);
    // setup the arrays that will store the interpolated and exact values at the Yee cell centres
    double ***Ex_exact = allocate3dmemory(Nx, Ny, Nz), ***Ex_interp = allocate3dmemory(Nx, Ny, Nz),
           ***Ey_exact = allocate3dmemory(Nx, Ny, Nz), ***Ey_interp = allocate3dmemory(Nx, Ny, Nz),
           ***Ez_exact = allocate3dmemory(Nx, Ny, Nz), ***Ez_interp = allocate3dmemory(Nx, Ny, Nz);

    // compute the exact field and the "split field" components
    // the interpolation functions are expecting split fields, but we can bypass this by making one split field component equal to _the entire field_ value, and the other zero
    for(int ii=0; ii<Nx; ii++) {
        for(int jj=0; jj<Ny; jj++) {
            for(int kk=0; kk<Nz; kk++) {
                // x components EX_{ii,jj,kk} := Ex[2*ii+1, 2*jj, 2*kk]
                Exy[kk][jj][ii] = Ex(x[2 * ii + 1], y[2 * jj], z[2 * kk]);
                Exz[kk][jj][ii] = 0.0;
                // y components EY_{ii,jj,kk} := Ey[2*ii, 2*jj+1, 2*kk]
                Eyx[kk][jj][ii] = Ey(x[2 * ii], y[2 * jj + 1], z[2 * kk]);
                Eyz[kk][jj][ii] = 0.0;
                // z components EZ_{ii,jj,kk} := Ez[2*ii, 2*jj, 2*kk+1]
                Ezx[kk][jj][ii] = Ez(x[2 * ii], y[2 * jj], z[2 * kk + 1]);
                Ezy[kk][jj][ii] = 0.0;

                // exact field components
                Ex_exact[kk][jj][ii] = Ex(x[2 * ii], y[2 * jj], z[2 * kk]);
                Ey_exact[kk][jj][ii] = Ey(x[2 * ii], y[2 * jj], z[2 * kk]);
                Ez_exact[kk][jj][ii] = Ez(x[2 * ii], y[2 * jj], z[2 * kk]);
            }
        }
    }

    // now we try to interpolate
    // for the sake of argument we assume a PML of 0 cells
    // so we need to interpolate to the centre of cells 1 through N{x,y,z}-1
    // note we still can't interpolate to the centre of cell 0 (until BAND_LIMITED_CELL_ZERO is implimented)
    for(int ii=1; ii<Nx; ii++) {
        for(int jj=1; jj<Ny; jj++) {
            for(int kk=1; kk<Nz; kk++) {
                // interpolate to the centre of cell (ii,jj,kk)
                double *Ex_val, *Ey_val, *Ez_val;
                interpolateTimeDomainEField(Exy, Exz, Eyx, Eyz, Ezx, Ezy,
                                            ii, jj, kk, Nx, Ny, Nz,
                                            &Ex_interp[kk][jj][ii], &Ey_interp[kk][jj][ii], &Ez_interp[kk][jj][ii]);
            }
        }
    }

    // can now deallocate our "fake" split field arrays, and "input data" arrays
    delete Exy, Exz, Eyx, Eyz, Ezx, Ezy;

    // values to hold the maximum errors
    double Ex_max_error = 0., Ey_max_error = 0., Ez_max_error = 0.;
    // now we compare the values across the exact and interp arrays,
    // across all indices from 1 to N{x,y,z}.
    for (int ii = 1; ii < Nx; ii++) {
        for (int jj = 1; jj < Ny; jj++) {
            for (int kk = 1; kk < Nz; kk++) {
                double curr_err_Ex = abs(Ex_exact[kk][jj][ii] - Ex_interp[kk][jj][ii]);
                if (curr_err_Ex > Ex_max_error){ Ex_max_error = curr_err_Ex;}

                double curr_err_Ey = abs(Ey_exact[kk][jj][ii] - Ey_interp[kk][jj][ii]);
                if (curr_err_Ey > Ey_max_error) { Ey_max_error = curr_err_Ey;}

                double curr_err_Ez = abs(Ez_exact[kk][jj][ii] - Ez_interp[kk][jj][ii]);
                if (curr_err_Ez > Ez_max_error) { Ez_max_error = curr_err_Ez;}
            }
        }
    }

    // delete large arrays to save memory again
    delete Ex_exact, Ex_interp, Ey_exact, Ey_interp, Ez_exact, Ez_interp;

    // check max errors do not exceed tolerances
    CHECK(Ex_max_error < Ex_tol);
    CHECK(Ey_max_error < Ey_tol);
    CHECK(Ez_max_error < Ez_tol);
}

/**
 * @brief Test the interpolation of the H-field components to the centre of the Yee cells
 * 
 */
TEST_CASE("H-field interpolation check") {

    // error tolerance
    // this needs to be set based on some kind of reference value.
    // Implementing this scheme in MATLAB (for the fields above) yields an error of
    // Hx : 6.35268084e-02 (using y-component 1st)
    // Hy : 9.99873462e-01 (using z-component 1st)
    // Hz : 9.99545033e-01 (using y-component 1st)
    double Hx_tol = 6.35268084e-02, Hy_tol = 9.99873462e-01, Hz_tol = 9.99545033e-01;

    // fake domain setup
    int Nx = 100, Ny = 50, Nz = 25;
    double Dx = 1. / (double)Nx, Dy = 1. / (double)Ny, Dz = 1. / (double)Nz;

    // setup the gridpoints
    double x[2 * Nx + 1], y[2 * Ny + 1], z[2 * Nz + 1];
    for (int i = 0; i < 2 * Nx + 1; i++)
    {
        x[i] = i * Dx / 2.;
    }
    for (int j = 0; j < 2 * Ny + 1; j++)
    {
        y[j] = j * Dy / 2.;
    }
    for (int k = 0; k < 2 * Nz + 1; k++)
    {
        z[k] = k * Dz / 2.;
    }

    // setup the "split" E-field components
    double ***Hxy = allocate3dmemory(Nx, Ny, Nz), ***Hxz = allocate3dmemory(Nx, Ny, Nz),
           ***Hyx = allocate3dmemory(Nx, Ny, Nz), ***Hyz = allocate3dmemory(Nx, Ny, Nz),
           ***Hzx = allocate3dmemory(Nx, Ny, Nz), ***Hzy = allocate3dmemory(Nx, Ny, Nz);
    // setup the arrays that will store the interpolated and exact values at the Yee cell centres
    double ***Hx_exact = allocate3dmemory(Nx, Ny, Nz), ***Hx_interp = allocate3dmemory(Nx, Ny, Nz),
           ***Hy_exact = allocate3dmemory(Nx, Ny, Nz), ***Hy_interp = allocate3dmemory(Nx, Ny, Nz),
           ***Hz_exact = allocate3dmemory(Nx, Ny, Nz), ***Hz_interp = allocate3dmemory(Nx, Ny, Nz);

    // compute the exact field and the "split field" components
    //Hx | (0.5, 0.0, 0.5)Hy | (0.5, 0.0, 0.5)Hz | (0.5, 0.5, 0.0)
    for (int ii = 0; ii < Nx; ii++)
    {
        for (int jj = 0; jj < Ny; jj++)
        {
            for (int kk = 0; kk < Nz; kk++)
            {
                // x components HX_{ii,jj,kk} := Hx[2*ii+1, 2*jj, 2*kk+1]
                Hxy[kk][jj][ii] = Hx(x[2 * ii + 1], y[2 * jj], z[2 * kk + 1]);
                Hxz[kk][jj][ii] = 0.;
                // y components HY_{ii,jj,kk} := Hy[2*ii+1, 2*jj, 2*kk+1]
                Hyx[kk][jj][ii] = Hy(x[2 * ii + 1], y[2 * jj], z[2 * kk + 1]);
                Hyz[kk][jj][ii] = 0.;
                // z components HZ_{ii,jj,kk} := Hz[2*ii+1, 2*jj+1, 2*kk]
                Hzx[kk][jj][ii] = Hz(x[2 * ii + 1], y[2 * jj + 1], z[2 * kk]);
                Hzy[kk][jj][ii] = 0.;

                // exact field components
                Hx_exact[kk][jj][ii] = Hx(x[2 * ii], y[2 * jj], z[2 * kk]);
                Hy_exact[kk][jj][ii] = Hy(x[2 * ii], y[2 * jj], z[2 * kk]);
                Hz_exact[kk][jj][ii] = Hz(x[2 * ii], y[2 * jj], z[2 * kk]);
            }
        }
    }

    // now we try to interpolate
    // for the sake of argument we assume a PML of 0 cells
    // so we need to interpolate to the centre of cells 1 through N{x,y,z}-1
    // note we still can't interpolate to the centre of cell 0 (until BAND_LIMITED_CELL_ZERO is implimented)
    for (int ii = 1; ii < Nx; ii++)
    {
        for (int jj = 1; jj < Ny; jj++)
        {
            for (int kk = 1; kk < Nz; kk++)
            {
                // interpolate to the centre of cell (ii,jj,kk)
                double *Hx_val, *Hy_val, *Hz_val;
                interpolateTimeDomainHField(Hxy, Hxz, Hyx, Hyz, Hzx, Hzy,
                                            ii, jj, kk, Nx, Ny, Nz,
                                            &Hx_interp[kk][jj][ii], &Hy_interp[kk][jj][ii], &Hz_interp[kk][jj][ii]);
            }
        }
    }

    // now we compare the values across the exact and interp arrays,
    // across all indices from 1 to N{x,y,z}.
    // cout << "ii, jj, kk \t | x-diff \t | y-diff \t | z-diff \n";
    for (int ii = 1; ii < Nx; ii++)
    {
        for (int jj = 1; jj < Ny; jj++)
        {
            for (int kk = 1; kk < Nz; kk++)
            {
                CHECK(abs(Hx_exact[kk][jj][ii] - Hx_interp[kk][jj][ii]) < Hx_tol);
                CHECK(abs(Hy_exact[kk][jj][ii] - Hy_interp[kk][jj][ii]) < Hy_tol);
                CHECK(abs(Hz_exact[kk][jj][ii] - Hz_interp[kk][jj][ii]) < Hz_tol);
                // cout << ii; cout << ", "; cout << jj; cout << ", "; cout << kk; cout << " \t | ";
                // cout << abs(Hx_exact[kk][jj][ii] - Hx_interp[kk][jj][ii]); cout << "\t | ";
                // cout << abs(Hy_exact[kk][jj][ii] - Hy_interp[kk][jj][ii]); cout << "\t | ";
                // cout << abs(Hz_exact[kk][jj][ii] - Hz_interp[kk][jj][ii]); cout << "\n";
            }
        }
    }
}