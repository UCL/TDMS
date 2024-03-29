#include "field.h"

#include "globals.h"
#include "interpolation_methods.h"

using namespace std;
using namespace tdms_math_constants;

double ElectricField::phase(int n, double omega, double dt) {
  return omega * ((double) n + 1) * dt;
}

void ElectricField::interpolate_transverse_electric_components(
        CellCoordinate cell, complex<double> *x_at_centre,
        complex<double> *y_at_centre, complex<double> *z_at_centre) {
  *x_at_centre = interpolate_to_centre_of(AxialDirection::X, cell);
  *y_at_centre = interpolate_to_centre_of(AxialDirection::Y, cell);
  *z_at_centre = complex<double>(0., 0.);
}
void ElectricField::interpolate_transverse_magnetic_components(
        CellCoordinate cell, complex<double> *x_at_centre,
        complex<double> *y_at_centre, complex<double> *z_at_centre) {
  *x_at_centre = complex<double>(0., 0.);
  *y_at_centre = complex<double>(0., 0.);
  *z_at_centre = interpolate_to_centre_of(AxialDirection::Z, cell);
}

complex<double> ElectricField::interpolate_to_centre_of(AxialDirection d,
                                                        CellCoordinate cell) {
  const InterpolationScheme *scheme;
  int i = cell.i, j = cell.j, k = cell.k;
  // prepare input data - if using a cubic scheme we have reserved more memory
  // than necessary but nevermind
  complex<double> interp_data[8];

  switch (d) {
    case X:
      // determine the interpolation scheme to use
      scheme = &(best_scheme(tot.i, i, interpolation_method));
      // now fill the interpolation data
      // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell
      // that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff;
           ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] =
                real.x[k][j][i - scheme->number_of_datapoints_to_left + ind] +
                IMAGINARY_UNIT *
                        imag.x[k][j]
                              [i - scheme->number_of_datapoints_to_left + ind];
      }
      break;
    case Y:
      // if we are in a 2D simulation, we just return the field value at cell
      // (i, 0, k) since there is no y-dimension to interpolate in.
      if (tot.j <= 1) {
        return complex<double>(real.y[k][0][i], imag.y[k][0][i]);
      } else {// 3D simulation, interpolation is as normal
        // determine the interpolation scheme to use
        scheme = &(best_scheme(tot.j, j, interpolation_method));

        // now fill the interpolation data
        // j - scheme.number_of_datapoints_to_left is the index of the Yee cell
        // that plays the role of v0 in the interpolation
        for (int ind = scheme->first_nonzero_coeff;
             ind <= scheme->last_nonzero_coeff; ind++) {
          interp_data[ind] =
                  real.y[k][j - scheme->number_of_datapoints_to_left + ind][i] +
                  IMAGINARY_UNIT *
                          imag.y[k][j - scheme->number_of_datapoints_to_left +
                                    ind][i];
        }
      }
      break;
    case Z:
      // determine the interpolation scheme to use
      scheme = &(best_scheme(tot.k, k, interpolation_method));

      // now fill the interpolation data
      // k - scheme.number_of_datapoints_to_left is the index of the Yee cell
      // that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff;
           ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] =
                real.z[k - scheme->number_of_datapoints_to_left + ind][j][i] +
                IMAGINARY_UNIT *
                        imag.z[k - scheme->number_of_datapoints_to_left + ind]
                              [j][i];
      }
      break;
    default:
      throw runtime_error(
              "Invalid axial direction selected for interpolation!\n");
      break;
  }
  // now run the interpolation scheme and place the result into the output
  return scheme->interpolate(interp_data);
}

double ElectricSplitField::interpolate_to_centre_of(AxialDirection d,
                                                    CellCoordinate cell) {
  const InterpolationScheme *scheme;
  int i = cell.i, j = cell.j, k = cell.k;
  // prepare input data - if using a cubic scheme we have reserved more memory
  // than necessary but nevermind
  double interp_data[8];

  switch (d) {
    case X:
      scheme = &(best_scheme(tot.i, i, interpolation_method));
      // now fill the interpolation data
      // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell
      // that plays the role of v0 in the interpolation
      for (int ind = scheme->first_nonzero_coeff;
           ind <= scheme->last_nonzero_coeff; ind++) {
        interp_data[ind] =
                xy(i - scheme->number_of_datapoints_to_left + ind, j, k) +
                xz(i - scheme->number_of_datapoints_to_left + ind, j, k);
      }
      // now run the interpolation scheme and place the result into the output
      return scheme->interpolate(interp_data);
      break;
    case Y:
      // if we are in a 2D simulation, we just return the field value at cell
      // (i, 0, k) since there is no y-dimension to interpolate in.
      if (tot.j <= 1) {
        return yx(i, 0, k) + yz(i, 0, k);
      } else {// 3D simulation, interpolation is as normal
        scheme = &(best_scheme(tot.j, j, interpolation_method));
        // now fill the interpolation data
        // j - scheme.number_of_datapoints_to_left is the index of the Yee cell
        // that plays the role of v0 in the interpolation
        for (int ind = scheme->first_nonzero_coeff;
             ind <= scheme->last_nonzero_coeff; ind++) {
          interp_data[ind] =
                  yx(i, j - scheme->number_of_datapoints_to_left + ind, k) +
                  yz(i, j - scheme->number_of_datapoints_to_left + ind, k);
        }
        // now run the interpolation scheme and place the result into the output
        return scheme->interpolate(interp_data);
        break;
        case Z:
          scheme = &(best_scheme(tot.k, k, interpolation_method));
          // now fill the interpolation data
          // k - scheme.number_of_datapoints_to_left is the index of the Yee
          // cell that plays the role of v0 in the interpolation
          for (int ind = scheme->first_nonzero_coeff;
               ind <= scheme->last_nonzero_coeff; ind++) {
            interp_data[ind] =
                    zx(i, j, k - scheme->number_of_datapoints_to_left + ind) +
                    zy(i, j, k - scheme->number_of_datapoints_to_left + ind);
          }
          // now run the interpolation scheme and place the result into the
          // output
          return scheme->interpolate(interp_data);
      }
      break;
    default:
      throw runtime_error(
              "Invalid axial direction selected for interpolation!\n");
      break;
  }
}
