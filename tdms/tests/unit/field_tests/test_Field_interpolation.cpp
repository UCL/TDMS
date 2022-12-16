/**
 * @file test_field_interpolation.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests interpolation of E- and H-fields and compares the errors against MATLAB benchmarks
 */
#include "field.h"

#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "cell_coordinate.h"
#include "globals.h"
#include "unit_test_utils.h"

using namespace tdms_math_constants;

using tdms_tests::TOLERANCE;
using tdms_tests::euclidean;
using tdms_tests::is_close_or_better;

/* Overview of Tests

We will test the performance of BLi at interpolating a known field to the centre of each Yee cell, for both the E- and H- fields.
The benchmark for success will be superior or equivalent error (Frobenius and 1D-dimension-wise norm) to that produced by MATLAB performing the same functionality.
Frobenious error is the Frobenius norm of the 3D matrix whose (i,j,k)-th entry is the error in the interpolated value at Yee cell centre (i,j,k).
The slice error, for a fixed j,k, is the norm-error of the 1D array of interpolated values along the axis (:,j,k).
The maximum of this is then the max-slice error: keeping track of this ensures us that the behaviour of BLi is consistent, and does not dramaically over-compensate in some areas and under-compensate in others.

All tests will be performed with cell sizes Dx = 0.25, Dy = 0.1, Dz = 0.05, over the range [-2,2].
*/

// functional form for the {E,H}-field components
inline double field_component(double t) {
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
  double Ex_fro_tol = 2.4535659128911650e-02, Ex_ms_tol = 1.0294466335592440e-03;
  double Ey_fro_tol = 2.5021563893394754e-02, Ey_ms_tol = 1.6590813667643383e-03;
  double Ez_fro_tol = 2.5181324617226587e-02, Ez_ms_tol = 1.7507103927894884e-03;

  // fake domain setup
  double x_lower = -2., y_lower = -2., z_lower = -2.;
  double extent_x = 4., extent_y = 4., extent_z = 4.;
  double cellDims[3] = {0.25, 0.1, 0.05};
  // The number of cells in each direction is then 16 = 4/0.25, 40 = 4/0.1, 80 = 4/0.05.
  // Note that due to the possibility that Nx/cellDims[0] computing something that is not quite an integer, we need to use round() to get an int safely
  int Nx = round(extent_x / cellDims[0]), Ny = round(extent_y / cellDims[1]),
      Nz = round(extent_z / cellDims[2]);

  // setup the "split" E-field components
  ElectricSplitField E_split(Nx-1, Ny-1, Nz-1);
  E_split.allocate(); // alocates Nx, Ny, Nz memory space here
  E_split.I_tot++; E_split.J_tot++; E_split.K_tot++; // correct the "number of datapoints" variable for these fields
  // setup for non-split field components
  ElectricField E(Nx, Ny, Nz);
  E.allocate();

  /* Compute the exact field and the "split field" components
  Ex[k][j][i] is the value of the field at position (x_lower,y_lower,z_lower) + (i,j+0.5,k+0.5)*cellDims
  Ey[k][j][i] is the value of the field at position (x_lower,y_lower,z_lower) + (i+0.5,j,k+0.5)*cellDims
  Ez[k][j][i] is the value of the field at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k)*cellDims
  */
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        double x_comp_value =
                field_component(y_lower + ((double) jj + 0.5) * cellDims[1]);// Ex depends on y
        double y_comp_value =
                field_component(z_lower + ((double) kk + 0.5) * cellDims[2]);// Ey depends on z
        double z_comp_value =
                field_component(x_lower + ((double) ii + 0.5) * cellDims[0]);// Ez depends on x
        // assign to "freq domain" ElectricField
        E.real.x[kk][jj][ii] = x_comp_value;
        E.imag.x[kk][jj][ii] = 0.;
        E.real.y[kk][jj][ii] = y_comp_value;
        E.imag.y[kk][jj][ii] = 0.;
        E.real.z[kk][jj][ii] = z_comp_value;
        E.imag.z[kk][jj][ii] = 0.;
        // assign to "time domain" ElectricSplitField - use weighting that sums to 1 to check addition is behaving as planned
        E_split.xy[kk][jj][ii] = x_comp_value;
        E_split.xz[kk][jj][ii] = 0.;
        E_split.yx[kk][jj][ii] = y_comp_value * .5;
        E_split.yz[kk][jj][ii] = y_comp_value * .5;
        E_split.zx[kk][jj][ii] = z_comp_value * .25;
        E_split.zy[kk][jj][ii] = z_comp_value * .75;
      }
    }
  }

  /* In each axis, we will now interpolate to positions 1 through N{x,y,z}-1. Recall this means to the midpoint of the datapoints indexed by 0 and 1, then 1 and 2, ..., N{x,y,z}-3 and N{x,y,z}-2, and finally N{x,y,z}-2 and N{x,y,z}-1.
  Ex_exact[k][j][i] is the field component at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k+0.5)*cellDims
  Ey_exact[k][j][i] is the field component at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k+0.5)*cellDims
  Ez_exact[k][j][i] is the field component at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k+0.5)*cellDims
  */
  Tensor3D<double> Ex_error, Ex_split_error, Ey_error, Ey_split_error, Ez_error, Ez_split_error;
  Ex_error.allocate(Nz, Ny, Nx-1);
  Ex_split_error.allocate(Nz, Ny, Nx-1);
  Ey_error.allocate(Nz, Ny-1, Nx);
  Ey_split_error.allocate(Nz, Ny-1, Nx);
  Ez_error.allocate(Nz-1, Ny, Nx);
  Ez_split_error.allocate(Nz-1, Ny, Nx);
  // now interpolate
  // note that we aren't interpolating to position 0 (before 1st point) or N{x,y,z} (after last point)
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        // coordinates of the cell that we are interested in
        double x_eval_position = y_lower + ((double) jj + 0.5) * cellDims[1];
        double y_eval_position = z_lower + ((double) kk + 0.5) * cellDims[2];
        double z_eval_position = x_lower + ((double) ii + 0.5) * cellDims[0];
        // current cell index
        CellCoordinate current_cell(ii, jj, kk);

        // Ex interpolation
        if (ii!=0) {
          double Ex_exact = field_component(x_eval_position);// Ex depends on y
          double Ex_split_interp =
                  E_split.interpolate_to_centre_of(AxialDirection::X, current_cell);
          double Ex_interp = E.interpolate_to_centre_of(AxialDirection::X, current_cell).real();
          Ex_error[kk][jj][ii-1] = Ex_interp - Ex_exact;
          Ex_split_error[kk][jj][ii-1] = Ex_split_interp - Ex_exact;
        }

        // Ey interpolation
        if (jj!=0) {
          double Ey_exact = field_component(y_eval_position);// Ey depends on z
          double Ey_split_interp =
                  E_split.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          double Ey_interp = E.interpolate_to_centre_of(AxialDirection::Y, current_cell).real();
          Ey_error[kk][jj-1][ii] = Ey_interp - Ey_exact;
          Ey_split_error[kk][jj-1][ii] = Ey_split_interp - Ey_exact;
        }

        // Ez interpolation
        if (kk!=0) {
          double Ez_exact = field_component(z_eval_position);// Ez depends on x
          double Ez_split_interp =
                  E_split.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          double Ez_interp = E.interpolate_to_centre_of(AxialDirection::Z, current_cell).real();
          Ez_error[kk-1][jj][ii] = Ez_interp - Ez_exact;
          Ez_split_error[kk-1][jj][ii] = Ez_split_interp - Ez_exact;
        }
      }
    }
  }
  // compute error-matrix Frobenius norms
  double Ex_fro_err = Ex_error.frobenius(), Ey_fro_err = Ey_error.frobenius(),
         Ez_fro_err = Ez_error.frobenius(), Ex_split_fro_err = Ex_split_error.frobenius(),
         Ey_split_fro_err = Ey_split_error.frobenius(),
         Ez_split_fro_err = Ez_split_error.frobenius();

  // compute max-slice errors
  double Ex_ms_err = 0., Ey_ms_err = 0., Ez_ms_err = 0., Ex_split_ms_err = 0., Ey_split_ms_err = 0.,
         Ez_split_ms_err = 0.;
  // Ex-slices
  for (int jj = 0; jj < Ny; jj++) {
    for (int kk = 0; kk < Nz; kk++) {
      // "slices" might not constitute sequential memory
      // as such, make a new array for safety
      double jk_errors[Nx - 1], jk_split_errors[Nx - 1];
      for (int ii = 0; ii < Nx - 1; ii++) {
        jk_errors[ii] = Ex_error[kk][jj][ii];
        jk_split_errors[ii] = Ex_split_error[kk][jj][ii];
      }
      // compute norm-error of this slice
      double jk_slice_error = euclidean(jk_errors, Nx - 1),
             jk_split_slice_error = euclidean(jk_split_errors, Nx - 1);
      // if this exceeds the current recorded maximum error, record this
      if (jk_slice_error > Ex_ms_err) { Ex_ms_err = jk_slice_error; }
      if (jk_split_slice_error > Ex_split_ms_err) { Ex_split_ms_err = jk_split_slice_error; }
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
      if (ik_slice_error > Ey_ms_err) { Ey_ms_err = ik_slice_error; }
      if (ik_split_slice_error > Ey_split_ms_err) { Ey_split_ms_err = ik_split_slice_error; }
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
      if (ij_slice_error > Ez_ms_err) { Ez_ms_err = ij_slice_error; }
      if (ij_split_slice_error > Ez_split_ms_err) { Ez_split_ms_err = ij_split_slice_error; }
    }
  }

  // check Frobenius errors are acceptable
  CHECK(is_close_or_better(Ex_fro_err, Ex_fro_tol));
  CHECK(is_close_or_better(Ey_fro_err, Ey_fro_tol));
  CHECK(is_close_or_better(Ez_fro_err, Ez_fro_tol));
  CHECK(is_close_or_better(Ex_split_fro_err, Ex_fro_tol));
  CHECK(is_close_or_better(Ey_split_fro_err, Ey_fro_tol));
  CHECK(is_close_or_better(Ez_split_fro_err, Ez_fro_tol));

  // check max-slice errors are acceptable
  CHECK(is_close_or_better(Ex_ms_err, Ex_ms_tol));
  CHECK(is_close_or_better(Ey_ms_err, Ey_ms_tol));
  CHECK(is_close_or_better(Ez_ms_err, Ez_ms_tol));
  CHECK(is_close_or_better(Ex_split_ms_err, Ex_ms_tol));
  CHECK(is_close_or_better(Ey_split_ms_err, Ey_ms_tol));
  CHECK(is_close_or_better(Ez_split_ms_err, Ez_ms_tol));

  // print information to the debugger/log
  SPDLOG_INFO(" Component | Frobenius err. : (  benchmark   ) | Max-slice err. : (  benchmark   )");
  SPDLOG_INFO("    x      | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ex_fro_err, Ex_fro_tol,
              Ex_ms_err, Ex_ms_tol);
  SPDLOG_INFO("    y      | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ey_fro_err, Ey_fro_tol,
              Ey_ms_err, Ey_ms_tol);
  SPDLOG_INFO("    z      | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ez_fro_err, Ez_fro_tol,
              Ez_ms_err, Ez_ms_tol);
  SPDLOG_INFO(" [Split] Component | Frobenius err. : (  benchmark   ) | Max-slice err. : (  "
              "benchmark   )");
  SPDLOG_INFO("        x          | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ex_split_fro_err,
              Ex_fro_tol, Ex_split_ms_err, Ex_ms_tol);
  SPDLOG_INFO("        y          | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ey_split_fro_err,
              Ey_fro_tol, Ey_split_ms_err, Ey_ms_tol);
  SPDLOG_INFO("        z          | {0:.8e} : ({1:.8e}) | {2:.8e} : ({3:.8e})", Ez_split_fro_err,
              Ez_fro_tol, Ez_split_ms_err, Ez_ms_tol);
}

/**
 * @brief Test the interpolation of the H-field components to the centre of the Yee cells
 *
 * Each component of the H-field will take the form
 * H_t(tt) = sin(2\pi tt) * exp(-tt^2)
 * The decision to make this function wholey imaginary is just to test whether the imaginary part of the field is correctly picked up and worked with.
 *
 * We only test Fro-norm error metrics, since interpolation must occur along two axes for each component.
 */
TEST_CASE("H-field interpolation check") {
  SPDLOG_INFO("===== Testing H-field BLi =====");
  // error tolerance, based on MATLAB performance
  // script: benchmark_test_field_interpolation_H.m
  double Hx_fro_tol = 6.1860732269207769e-02;
  double Hy_fro_tol = 1.2622903101481411e-01;
  double Hz_fro_tol = 7.3136417157073502e-02;

  // fake domain setup
  double x_lower = -2., y_lower = -2., z_lower = -2.;
  double extent_x = 4., extent_y = 4., extent_z = 4.;
  double cellDims[3] = {0.1, 0.05, 0.025};
  // The number of cells in each direction is then 40, 80,  160 respectively
  // Note that due to the possibility that Nx/cellDims[0] computing something that is not quite an integer, we need to use round() to get an int safely
  int Nx = round(extent_x / cellDims[0]), Ny = round(extent_y / cellDims[1]),
      Nz = round(extent_z / cellDims[2]);

  // setup the "split" H-field components
  MagneticSplitField H_split(Nx - 1, Ny - 1, Nz - 1);
  H_split.allocate();// alocates Nx, Ny, Nz memory space here
  // setup the non-split field components
  H_split.I_tot++;
  H_split.J_tot++;
  H_split.K_tot++;// correct the "number of datapoints" variable for these fields
  MagneticField H(Nx, Ny, Nz);
  H.allocate();

  /* Compute the exact field and the "split field" components
  Hx_exact[k-1][j-1][i] is the value of the field at position (x_lower,y_lower,z_lower) + (i+0.5,j,k)*cellDims
  Hy_exact[k-1][j][i-1] is the value of the field at position (x_lower,y_lower,z_lower) + (i,j+0.5,k)*cellDims
  Hz_exact[k][j-1][i-1] is the value of the field at position (x_lower,y_lower,z_lower) + (i,j,k+0.5)*cellDims
  */
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        double x_comp_value =
                field_component(y_lower + ((double) jj + 1.) * cellDims[1]);// Hx depends on y
        double y_comp_value =
                field_component(z_lower + ((double) kk + 1.) * cellDims[2]);// Hy depends on z
        double z_comp_value =
                field_component(x_lower + ((double) ii + 1.) * cellDims[0]);// Hz depends on x
        // assign to "freq domain" ElectricField
        H.imag.x[kk][jj][ii] = x_comp_value;
        H.real.x[kk][jj][ii] = 0.;
        H.imag.y[kk][jj][ii] = y_comp_value;
        H.real.y[kk][jj][ii] = 0.;
        H.imag.z[kk][jj][ii] = z_comp_value;
        H.real.z[kk][jj][ii] = 0.;
        // assign to "time domain" ElectricSplitField - use weighting that sums to 1 to check addition is behaving as planned
        H_split.xy[kk][jj][ii] = x_comp_value;
        H_split.xz[kk][jj][ii] = 0.;
        H_split.yx[kk][jj][ii] = y_comp_value * .25;
        H_split.yz[kk][jj][ii] = y_comp_value * .75;
        H_split.zx[kk][jj][ii] = z_comp_value * .125;
        H_split.zy[kk][jj][ii] = z_comp_value * .875;
      }
    }
  }

  /* In each axis, we will now interpolate to positions 1 through N{x,y,z}-1. Recall this means to the midpoint of the datapoints indexed by 0 and 1, then 1 and 2, ..., N{x,y,z}-3 and N{x,y,z}-2, and finally N{x,y,z}-2 and N{x,y,z}-1.
  Ex_exact[k][j][i] is the field component at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k+0.5)*cellDims
  Ey_exact[k][j][i] is the field component at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k+0.5)*cellDims
  Ez_exact[k][j][i] is the field component at position (x_lower,y_lower,z_lower) + (i+0.5,j+0.5,k+0.5)*cellDims
  */
  Tensor3D<double> Hx_error, Hx_split_error, Hy_error, Hy_split_error, Hz_error, Hz_split_error;
  Hx_error.allocate(Nz - 1, Ny - 1, Nx);
  Hx_split_error.allocate(Nz - 1, Ny - 1, Nx);
  Hy_error.allocate(Nz - 1, Ny, Nx - 1);
  Hy_split_error.allocate(Nz - 1, Ny, Nx - 1);
  Hz_error.allocate(Nz, Ny - 1, Nx - 1);
  Hz_split_error.allocate(Nz, Ny - 1, Nx - 1);

  // now interpolate
  // note that we aren't interpolating to position 0 (before 1st point) or N{x,y,z} (after last point)
  for (int ii = 0; ii < Nx; ii++) {
    for (int jj = 0; jj < Ny; jj++) {
      for (int kk = 0; kk < Nz; kk++) {
        // coordinates of the cell that we are interested in
        double x_eval_position = y_lower + ((double) jj + 0.5) * cellDims[1];
        double y_eval_position = z_lower + ((double) kk + 0.5) * cellDims[2];
        double z_eval_position = x_lower + ((double) ii + 0.5) * cellDims[0];
        // current cell index
        CellCoordinate current_cell(ii, jj, kk);

        // Hx interpolation
        if (jj != 0 && kk != 0) {
          double Hx_exact = field_component(x_eval_position);// Hx depends on y
          double Hx_split_interp =
                  H_split.interpolate_to_centre_of(AxialDirection::X, current_cell);
          double Hx_interp = H.interpolate_to_centre_of(AxialDirection::X, current_cell).imag();
          Hx_error[kk-1][jj-1][ii] = Hx_interp - Hx_exact;
          Hx_split_error[kk-1][jj-1][ii] = Hx_split_interp - Hx_exact;
        }

        // Hy interpolation
        if (ii != 0 && kk != 0) {
          double Hy_exact = field_component(y_eval_position);// Hy depends on z
          double Hy_split_interp =
                  H_split.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          double Hy_interp = H.interpolate_to_centre_of(AxialDirection::Y, current_cell).imag();
          Hy_error[kk-1][jj][ii-1] = Hy_interp - Hy_exact;
          Hy_split_error[kk-1][jj][ii-1] = Hy_split_interp - Hy_exact;
        }

        // Hz interpolation
        if (ii != 0 && jj != 0) {
          double Hz_exact = field_component(z_eval_position);// Hz depends on x
          double Hz_split_interp =
                  H_split.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          double Hz_interp = H.interpolate_to_centre_of(AxialDirection::Z, current_cell).imag();
          Hz_error[kk][jj-1][ii-1] = Hz_interp - Hz_exact;
          Hz_split_error[kk][jj-1][ii-1] = Hz_split_interp - Hz_exact;
        }
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
  CHECK(is_close_or_better(Hx_fro_err, Hx_fro_tol));
  CHECK(is_close_or_better(Hy_fro_err, Hy_fro_tol));
  CHECK(is_close_or_better(Hz_fro_err, Hz_fro_tol));
  CHECK(is_close_or_better(Hx_split_fro_err, Hx_fro_tol));
  CHECK(is_close_or_better(Hy_split_fro_err, Hy_fro_tol));
  CHECK(is_close_or_better(Hz_split_fro_err, Hz_fro_tol));

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
