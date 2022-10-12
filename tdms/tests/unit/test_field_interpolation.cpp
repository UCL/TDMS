/**
 * @file test_field_interpolation.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests interpolation of E- and H-fields and compares the errors against MATLAB benchmarks
 */
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <spdlog/spdlog.h>

#include "globals.h"
#include "field.h"

using namespace tdms_math_constants;

/* Overview of Tests

We will test the performance of BLi at interpolating a known field to the centre of each Yee cell, for both the E- and H- fields.
The benchmark for success will be superior or equivalent error (Frobenius and 1D-dimension-wise norm) to that produced by MATLAB performing the same functionality.
Frobenious error is the Frobenius norm of the 3D matrix whose (i,j,k)-th entry is the error in the interpolated value at Yee cell centre (i,j,k).
The slice error, for a fixed j,k, is the norm-error of the 1D array of interpolated values along the axis (:,j,k).
The maximum of this is then the max-slice error: keeping track of this ensures us that the behaviour of BLi is consistent, and does not dramaically over-compensate in some areas and under-compensate in others.

All tests will be performed with cell sizes Dx = 0.25, Dy = 0.1, Dz = 0.05, over the range [-2,2].
*/

// computes the Euclidean norm of a 1d-array (or pointer thereto)
inline double euclidean(double *v, int end, int start = 0) {
    double norm_val = 0.;
    for (int i = start; i < end; i++) {
        norm_val += v[i] * v[i];
    }
    return sqrt(norm_val);
}

// functional form for the E-field components
inline double Ecomponent(double t) {
    return sin(2. * DCPI * t) * exp(-t * t);
}
// functional form for the H-field components
inline double Hcomponent(double t) {
    return sin(2. * DCPI * t) * exp(-t * t);
}

/**
 * @brief Test the interpolation of the E-field components to the centre of the Yee cells
 *
 * Each component of the E-field will take the form
 * E_t(tt) = sin(2\pi tt) * exp(-tt^2)
 *
 * We test both the Fro- and slice-norm metrics, since interpolation only happens along one axis
 */
TEST_CASE("E-field interpolation check") {
    SPDLOG_INFO("===== Testing E-field BLi =====");
    // error tolerance, based on MATLAB performance
    // script: benchmark_test_field_interpolation_H.m
    double Ex_fro_tol = 1.8602247952705870e+00, Ex_ms_tol = 3.2884439181679728e-02;
    double Ey_fro_tol = 2.8009859776421876e-02, Ey_ms_tol = 7.8289938125395461e-04;
    double Ez_fro_tol = 1.0884557082700464e-02, Ez_ms_tol = 4.3024989629313981e-04;

    // additional tolerance to allow for floating-point rounding imprecisions, etc
    double acc_tol = 1e-12;

    // fake domain setup
    double x_lower = -2., y_lower = -2., z_lower = -2.;
    double extent_x = 4., extent_y = 4., extent_z = 4.;
    double cellDims[3] = {0.25, 0.1, 0.05};
    // The number of cells in each direction is then 16 = 4/0.25, 40 = 4/0.1, 80 = 4/0.05.
    // Note that due to the possibility that Nx/cellDims[0] computing something that is not quite an integer, we need to use round() to get an int safely
    int Nx = round(extent_x / cellDims[0]), Ny = round(extent_y / cellDims[1]), Nz = round(extent_z / cellDims[2]);
    SPDLOG_INFO("(Nx, Ny, Nz) = ({},{},{})", Nx, Ny, Nz);

    // setup the "split" E-field components
    ElectricSplitField E_split(Nx, Ny, Nz);
    E_split.allocate();
    // setup for non-split field components
    ElectricField E(Nx, Ny, Nz);
    E.allocate();
    // storage for pointwise errors (-1 since we don't interpolate to cell 0's centre)
    Tensor3D<double> Ex_error, Ex_split_error, Ey_error, Ey_split_error, Ez_error, Ez_split_error;
    Ex_error.allocate(Nz, Ny, Nx - 1); Ex_split_error.allocate(Nz, Ny, Nx - 1);
    Ey_error.allocate(Nz, Ny - 1, Nx); Ey_split_error.allocate(Nz, Ny - 1, Nx);
    Ez_error.allocate(Nz - 1, Ny, Nx); Ez_split_error.allocate(Nz - 1, Ny, Nx);

    // compute the exact field and the "split field" components
    // the interpolation functions are expecting split fields, but we can bypass this by making one split field component equal to _the entire field_ value, and the other zero
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                // these are the coordinates of Yee cell i,j,k's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // Initialise sample values that we will pass to the interpolation schemes:
                // in each case, set one component of the split field to be the "whole" field,
                // and the other to be 0.

                // E{x,y,z} offsets from cell centres are 0.5*D{x,y,z}
                double x_comp_value = Ecomponent(cell_centre[0] + 0.5 * cellDims[0]),
                       y_comp_value = Ecomponent(cell_centre[1] + 0.5 * cellDims[1]),
                       z_comp_value = Ecomponent(cell_centre[2] + 0.5 * cellDims[2]);
                // assign component values
                E.real.x[kk][jj][ii] = x_comp_value; E.imag.x[kk][jj][ii] = 0.;
                E.real.y[kk][jj][ii] = y_comp_value; E.imag.y[kk][jj][ii] = 0.;
                E.real.z[kk][jj][ii] = z_comp_value; E.imag.z[kk][jj][ii] = 0.;
                // split fields - use some wieghting that sums to one for the split cells
                E_split.xy[kk][jj][ii] = x_comp_value;
                E_split.xz[kk][jj][ii] = 0.;
                E_split.yx[kk][jj][ii] = y_comp_value * .5;
                E_split.yz[kk][jj][ii] = y_comp_value * .5;
                E_split.zx[kk][jj][ii] = z_comp_value * .25;
                E_split.zy[kk][jj][ii] = z_comp_value * .75;
            }
        }
    }

    // run interpolation functions

    // interpolate Ex
    for (int ii = 0; ii < Nx - 1; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
          for (int kk = 0; kk < Nz; kk++) {
            // we are interpolating to the centre of cell ii+1,jj,kk
            // these are the coordinates of the Yee cell's centre
            double cell_centre[3];
            cell_centre[0] = x_lower + ((double) ii + 1.5) * cellDims[0];
            cell_centre[1] = y_lower + ((double) jj + 0.5) * cellDims[1];
            cell_centre[2] = z_lower + ((double) kk + 0.5) * cellDims[2];

            // compute the true value of the field components at the centre of this Yee cell
            double Ex_exact = Ecomponent(cell_centre[0]);

            // interpolate to the centre of this cell
            double Ex_interp, Ex_split_interp;
            Ex_interp = E.interpolate_x_to_centre(ii + 1, jj, kk).real();
            Ex_split_interp = E_split.interpolate_x_to_centre(ii + 1, jj, kk);

            // compute the errors
            Ex_error[kk][jj][ii] = Ex_interp - Ex_exact;
            Ex_split_error[kk][jj][ii] = Ex_split_interp - Ex_exact;
          }
        }
    }
    // interpolate Ey
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny - 1; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                // we are interpolating to the centre of cell ii,jj+1,kk
                // these are the coordinates of the Yee cell's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 1.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Ey_exact = Ecomponent(cell_centre[1]);

                // interpolate to the centre of this cell
                double Ey_interp, Ey_split_interp;
                Ey_interp = E.interpolate_y_to_centre(ii, jj + 1, kk).real();
                Ey_split_interp = E_split.interpolate_y_to_centre(ii, jj + 1, kk);

                // compute the errors
                Ey_error[kk][jj][ii] = Ey_interp - Ey_exact;
                Ey_split_error[kk][jj][ii] = Ey_split_interp - Ey_exact;
            }
        }
    }
    // interpolate Ez
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            for (int kk = 0; kk < Nz - 1; kk++) {
                // we are interpolating to the centre of cell ii,jj,kk+1
                // these are the coordinates of the Yee cell's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 1.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Ez_exact = Ecomponent(cell_centre[2]);

                // interpolate to the centre of this cell
                double Ez_interp, Ez_split_interp;
                Ez_interp = E.interpolate_z_to_centre(ii, jj, kk + 1).real();
                Ez_split_interp = E_split.interpolate_z_to_centre(ii, jj, kk + 1);

                // compute the errors
                Ez_error[kk][jj][ii] = Ez_interp - Ez_exact;
                Ez_split_error[kk][jj][ii] = Ez_split_interp - Ez_exact;
            }
        }
    }

    // compute error-matrix Frobenius norms
    double Ex_fro_err = Ex_error.frobenius(),
           Ey_fro_err = Ey_error.frobenius(),
           Ez_fro_err = Ez_error.frobenius(),
           Ex_split_fro_err = Ex_split_error.frobenius(),
           Ey_split_fro_err = Ey_split_error.frobenius(),
           Ez_split_fro_err = Ez_split_error.frobenius();

    // compute max-slice errors
    double Ex_ms_err = 0., Ey_ms_err = 0., Ez_ms_err = 0.,
           Ex_split_ms_err = 0., Ey_split_ms_err = 0., Ez_split_ms_err = 0.;
    // Ex-slices
    for (int jj = 0; jj < Ny; jj++) {
        for (int kk = 0; kk < Nz; kk++) {
            // "slices" might not constitute sequential memory
            // as such, make a new array for safety
            double jk_errors[Nx - 1], jk_split_errors[Nx -1];
            for (int ii = 0; ii < Nx - 1; ii++) {
                jk_errors[ii] = Ex_error[kk][jj][ii];
                jk_split_errors[ii] = Ex_split_error[kk][jj][ii];
            }
            // compute norm-error of this slice
            double jk_slice_error = euclidean(jk_errors, Nx - 1),
                   jk_split_slice_error = euclidean(jk_split_errors, Nx - 1);
            // if this exceeds the current recorded maximum error, record this
            if (jk_slice_error > Ex_ms_err) {
                Ex_ms_err = jk_slice_error;
            }
            if (jk_split_slice_error > Ex_split_ms_err) {
                Ex_split_ms_err = jk_split_slice_error;
            }
        }
    }
    // Ey-slices
    for (int ii = 0; ii < Nx; ii++) {
        for (int kk = 0; kk < Nz; kk++) {
            double ik_errors[Ny - 1], ik_split_errors[Ny - 1];
            for (int jj = 0; jj < Ny - 1; jj++) {
                ik_errors[jj] = Ey_error[kk][jj][ii];
                ik_split_errors[jj] = Ey_split_error[kk][jj][ii];
            }
            double ik_slice_error = euclidean(ik_errors, Ny - 1),
                   ik_split_slice_error = euclidean(ik_split_errors, Ny - 1);
            if (ik_slice_error > Ey_ms_err) {
                Ey_ms_err = ik_slice_error;
            }
            if (ik_split_slice_error > Ey_split_ms_err) {
                Ey_split_ms_err = ik_split_slice_error;
            }
        }
    }
    // Ez-slices
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            double ij_errors[Nz - 1], ij_split_errors[Nz - 1];
            for (int kk = 0; kk < Nz - 1; kk++) {
                ij_errors[kk] = Ez_error[kk][jj][ii];
                ij_split_errors[kk] = Ez_split_error[kk][jj][ii];
            }
            double ij_slice_error = euclidean(ij_errors, Nz - 1),
                   ij_split_slice_error = euclidean(ij_split_errors, Nz - 1);
            if (ij_slice_error > Ez_ms_err) {
                Ez_ms_err = ij_slice_error;
            }
            if (ij_split_slice_error > Ez_split_ms_err) {
                Ez_split_ms_err = ij_split_slice_error;
            }
        }
    }

    // check Frobenius errors are acceptable
    CHECK(Ex_fro_err <= Ex_fro_tol + acc_tol);
    CHECK(Ey_fro_err <= Ey_fro_tol + acc_tol);
    CHECK(Ez_fro_err <= Ez_fro_tol + acc_tol);
    CHECK(Ex_split_fro_err <= Ex_fro_tol + acc_tol);
    CHECK(Ey_split_fro_err <= Ey_fro_tol + acc_tol);
    CHECK(Ez_split_fro_err <= Ez_fro_tol + acc_tol);

    // check max-slice errors are acceptable
    CHECK(Ex_ms_err <= Ex_ms_tol + acc_tol);
    CHECK(Ey_ms_err <= Ey_ms_tol + acc_tol);
    CHECK(Ez_ms_err <= Ez_ms_tol + acc_tol);
    CHECK(Ex_split_ms_err <= Ex_ms_tol + acc_tol);
    CHECK(Ey_split_ms_err <= Ey_ms_tol + acc_tol);
    CHECK(Ez_split_ms_err <= Ez_ms_tol + acc_tol);

    // print information to the debugger/log
    SPDLOG_INFO(" Component | Frobenius err. : (  benchmark   ) | Max-slice err. : (  benchmark   )");
    SPDLOG_INFO("    x      | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ex_fro_err, Ex_fro_tol, Ex_ms_err, Ex_ms_tol);
    SPDLOG_INFO("    y      | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ey_fro_err, Ey_fro_tol, Ey_ms_err, Ey_ms_tol);
    SPDLOG_INFO("    z      | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ez_fro_err, Ez_fro_tol, Ez_ms_err, Ez_ms_tol);
    SPDLOG_INFO(" [Split] Component | Frobenius err. : (  benchmark   ) | Max-slice err. : (  benchmark   )");
    SPDLOG_INFO("        x          | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ex_split_fro_err, Ex_fro_tol, Ex_split_ms_err, Ex_ms_tol);
    SPDLOG_INFO("        y          | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ey_split_fro_err, Ey_fro_tol, Ey_split_ms_err, Ey_ms_tol);
    SPDLOG_INFO("        z          | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ez_split_fro_err, Ez_fro_tol, Ez_split_ms_err, Ez_ms_tol);
}

/**
 * @brief Test the interpolation of the H-field components to the centre of the Yee cells
 *
 * Each component of the H-field will take the form
 * H_t(tt) = \mathrm{i} * ( sin(2\pi tt) * exp(-tt^2) ).
 * The decision to make this function wholey imaginary is just to test whether the imaginary part of the field is correctly picked up and worked with.
 *
 * We only test Fro-norm error metrics, since interpolation must occur along two axes for each component.
 */
TEST_CASE("H-field interpolation check") {
    SPDLOG_INFO("===== Testing H-field BLi =====");
    // error tolerance, based on MATLAB performance
    // script: benchmark_test_field_interpolation_H.m
    double Hx_fro_tol = 4.9722096274682832e-02;
    double Hy_fro_tol = 4.8756759347066955e-02;
    double Hz_fro_tol = 4.8447195218466632e-02;

    // additional tolerance to allow for floating-point rounding imprecisions, etc
    double acc_tol = 1e-12;

    // fake domain setup
    double x_lower = -2., y_lower = -2., z_lower = -2.;
    double extent_x = 4., extent_y = 4., extent_z = 4.;
    double cellDims[3] = {0.25, 0.1, 0.05};
    // The number of cells in each direction is then 16 = 4/0.25, 40 = 4/0.1, 80 = 4/0.05.
    // Note that due to the possibility that Nx/cellDims[0] computing something that is not quite an integer, we need to use round() to get an int safely
    int Nx = round(extent_x / cellDims[0]), Ny = round(extent_y / cellDims[1]), Nz = round(extent_z / cellDims[2]);
    SPDLOG_INFO("(Nx, Ny, Nz) = ({},{},{})", Nx, Ny, Nz);

    // setup the "split" H-field components
    MagneticSplitField H_split(Nx, Ny, Nz);
    H_split.allocate();
    // setup the non-split field components
    MagneticField H(Nx, Ny, Nz);
    H.allocate();
    // storage for pointwise errors (-1 since we don't interpolate to cell 0's centre)
    Tensor3D<double> Hx_error, Hx_split_error, Hy_error, Hy_split_error, Hz_error, Hz_split_error;
    Hx_error.allocate(Nz - 1, Ny - 1, Nx); Hx_split_error.allocate(Nz - 1, Ny - 1, Nx);
    Hy_error.allocate(Nz - 1, Ny, Nx - 1); Hy_split_error.allocate(Nz - 1, Ny, Nx - 1);
    Hz_error.allocate(Nz, Ny - 1, Nx - 1); Hz_split_error.allocate(Nz, Ny - 1, Nx - 1);

    // compute the exact field and the "split field" components
    // the interpolation functions are expecting split fields, but we can bypass this by making one split field component equal to _the entire field_ value, and the other zero
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                // these are the coordinates of Yee cell i,j,k's centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // Initialise sample values that we will pass to the interpolation schemes:
                // in each case, set one component of the split field to be the "whole" field,
                // and the other to be 0.

                // H{x,y,z} offsets from cell centres are 0.5 * D{!{x,y,z}} IE, 0.5 away from the centre in the two directions that are _not_ the field component
                double x_comp_value = Hcomponent(cell_centre[0]),
                       y_comp_value = Hcomponent(cell_centre[1]),
                       z_comp_value = Hcomponent(cell_centre[2]);
                // assign fields. For split field, use weightings that sum to unity
                H.imag.x[kk][jj][ii] = x_comp_value; H.real.x[kk][jj][ii] = 0.;
                H.imag.y[kk][jj][ii] = y_comp_value; H.real.y[kk][jj][ii] = 0.;
                H.imag.z[kk][jj][ii] = z_comp_value; H.real.z[kk][jj][ii] = 0.;
                H_split.xy[kk][jj][ii] = x_comp_value;
                H_split.xz[kk][jj][ii] = 0.;
                H_split.yx[kk][jj][ii] = .125 * y_comp_value;
                H_split.yz[kk][jj][ii] = .875 * y_comp_value;
                H_split.zx[kk][jj][ii] = .0625 * z_comp_value;
                H_split.zy[kk][jj][ii] = .9375 * z_comp_value;
            }
        }
    }

    // run interpolation functions

    // interpolate Hx
    for (int ii = 0; ii < Nx; ii++) {
        for (int jj = 0; jj < Ny - 1; jj++) {
            for (int kk = 0; kk < Nz - 1; kk++) {
                // we are interpolating to the centre of cell ii,jj+1,kk+1
                // these are the coordinates of Yee cell centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 0.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 1.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 1.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Hx_exact = Hcomponent(cell_centre[0]);

                // interpolate to the centre of this cell
                double Hx_interp, Hx_split_interp;
                Hx_interp = H.interpolate_x_to_centre(ii, jj + 1, kk + 1).imag();
                Hx_split_interp = H_split.interpolate_x_to_centre(ii, jj + 1, kk + 1);

                // compute the errors
                Hx_error[kk][jj][ii] = Hx_interp - Hx_exact;
                Hx_split_error[kk][jj][ii] = Hx_split_interp - Hx_exact;
            }
        }
    }
    // interpolate Hy
    for (int ii = 0; ii < Nx - 1; ii++) {
        for (int jj = 0; jj < Ny; jj++) {
            for (int kk = 0; kk < Nz - 1; kk++) {
                // we are interpolating to the centre of cell ii,jj+1,kk+1
                // these are the coordinates of Yee cell centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 1.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 0.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 1.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Hy_exact = Hcomponent(cell_centre[1]);

                // interpolate to the centre of this cell
                double Hy_interp, Hy_split_interp;
                Hy_interp = H.interpolate_y_to_centre(ii + 1, jj, kk + 1).imag();
                Hy_split_interp = H_split.interpolate_y_to_centre(ii + 1, jj, kk + 1);

                // compute the errors
                Hy_error[kk][jj][ii] = Hy_interp - Hy_exact;
                Hy_split_error[kk][jj][ii] = Hy_split_interp - Hy_exact;
            }
        }
    }
    // interpolate Hz
    for (int ii = 0; ii < Nx - 1; ii++) {
        for (int jj = 0; jj < Ny - 1; jj++) {
            for (int kk = 0; kk < Nz; kk++) {
                // we are interpolating to the centre of cell ii+1,jj+1,kk
                // these are the coordinates of Yee cell centre
                double cell_centre[3];
                cell_centre[0] = x_lower + ((double)ii + 1.5) * cellDims[0];
                cell_centre[1] = y_lower + ((double)jj + 1.5) * cellDims[1];
                cell_centre[2] = z_lower + ((double)kk + 0.5) * cellDims[2];

                // compute the true value of the field components at the centre of this Yee cell
                double Hz_exact = Hcomponent(cell_centre[2]);

                // interpolate to the centre of this cell
                double Hz_interp, Hz_split_interp;
                Hz_interp = H.interpolate_z_to_centre(ii + 1, jj + 1, kk).imag();
                Hz_split_interp = H_split.interpolate_z_to_centre(ii + 1, jj + 1, kk);

                // compute the errors
                Hz_error[kk][jj][ii] = Hz_interp - Hz_exact;
                Hz_split_error[kk][jj][ii] = Hz_split_interp - Hz_exact;
            }
        }
    }

    // compute Frobenius norms
    double Hx_fro_err = Hx_error.frobenius(),
           Hy_fro_err = Hy_error.frobenius(),
           Hz_fro_err = Hz_error.frobenius(),
           Hx_split_fro_err = Hx_split_error.frobenius(),
           Hy_split_fro_err = Hy_split_error.frobenius(),
           Hz_split_fro_err = Hz_split_error.frobenius();

    // check Frobenius errors are acceptable
    CHECK(Hx_fro_err <= Hx_fro_tol + acc_tol);
    CHECK(Hy_fro_err <= Hy_fro_tol + acc_tol);
    CHECK(Hz_fro_err <= Hz_fro_tol + acc_tol);
    CHECK(Hx_split_fro_err <= Hx_fro_tol + acc_tol);
    CHECK(Hy_split_fro_err <= Hy_fro_tol + acc_tol);
    CHECK(Hz_split_fro_err <= Hz_fro_tol + acc_tol);

    // print information to the debugger/log
    SPDLOG_INFO(" Component | Frobenius err. : (  benchmark   )");
    SPDLOG_INFO("    x      | {0:.8e} : ({1:.8e})", Hx_fro_err, Hx_fro_tol);
    SPDLOG_INFO("    y      | {0:.8e} : ({1:.8e})", Hy_fro_err, Hy_fro_tol);
    SPDLOG_INFO("    z      | {0:.8e} : ({1:.8e})", Hz_fro_err, Hz_fro_tol);
    SPDLOG_INFO(" [Split] Component | Frobenius err. : (  benchmark   )");
    SPDLOG_INFO("        x          | {0:.8e} : ({1:.8e})", Hx_split_fro_err, Hx_fro_tol);
    SPDLOG_INFO("        y          | {0:.8e} : ({1:.8e})", Hy_split_fro_err, Hy_fro_tol);
    SPDLOG_INFO("        z          | {0:.8e} : ({1:.8e})", Hz_split_fro_err, Hz_fro_tol);
}
