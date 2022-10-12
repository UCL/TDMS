#include "field.h"
#include "interpolation_methods.h"

using namespace std;


double MagneticField::phase(int n, double omega, double dt){
  return omega * ((double) n + 0.5) * dt;  // 0.5 added because it's known half a time step after E
}

/* 2D INTERPOLATION SCHEMES (FOR THE MAGNETIC FIELD IN 3D SIMULATIONS)

Unlike the E-field, the H-field components associated with Yee cell i,j,k are _not_ aligned with the centre of the Yee cells.
Instead, the position of the field components (relative to the Yee cell centre is):

Hx  | (0.0, 0.5, 0.5) .* (Dx, Dy, Dx)
Hy  | (0.5, 0.0, 0.5) .* (Dx, Dy, Dz)
Hz  | (0.5, 0.5, 0.0) .* (Dx, Dy, Dz)

where Dx, Dy, Dz are the dimensions of the Yee cell.
This requires us to interpolate twice to recover (any of the) field components at the centre of Yee cell i,j,k, when running a 3D simulation.

Henceforth, let {a,b,c} = {x,y,z} be some 1-to-1 assignment of the co-ordinate axes to the labels a,b,c.
The values Da, Db, Dc are the corresponding permutation of Dx, Dy, Dz.
Suppose that we wish to interpolate the Ha field component to the centre a particular Yee cell.

We will use the notation (a_i,b_j,c_k) for the Yee cell indices, although bear in mind that this does not reflect the order the indices appear in the code.
For example; if a = y, b = x, c = z, then the value of Ha at cell (a_i,b_j,c_k) is accessed via Ha[c_k][a_i][b_j], due to the interchanging of the x and y directions.
Similarly; we will write Ha[a_i, b_j, c_k] to refer to the value of Ha associated to cell (a_i,b_j,c_k).

Suppose now that we have selected cell a_i,b_j,c_k to interpolate to the centre of.
We must interpolate in the b-direction, then c-direction, or vice-versa.
For optimal accuracy we must determine the best interpolation scheme we can use in each direction at cell (a_i,b_j,c_k), and use the WORSE scheme second.

Let us assume WLOG that the b-direction interpolation scheme (b_scheme) is inferior to that of the c-direction (c_scheme).
We now need to interpolate Ha in the c-direction to obtain the value of Ha at the spatial positions (a_i, cell_b + Db, c_k) where
b_j - b_scheme.number_of_datapoints_to_left + b_scheme.first_nonzero_coeff <= cell_b <= b_j - b_scheme.number_of_datapoints_to_left + b_scheme.last_nonzero_coeff.

    Let cell_b be one particular index in this range.
    To apply c_scheme to obtain the value of Ha at (a_i, cell_b + Db, c_k), we require the values Ha[a_i, cell_b, cell_c] where
    c_k - c_scheme.number_of_datapoints_to_left + c_scheme.first_nonzero_coeff <= cell_c <= c_k - c_scheme.number_of_datapoints_to_left + c_scheme.last_nonzero_coeff.
    These values are fed to c_scheme.interpolate, which provides us with Ha at the spatial position (a_i, cell_b + Db, c_k).

Now with approximations of Ha at each (a_i, cell_b + Db, c_k), we can pass this information to b_scheme.interpolate to recover the value of Ha at the centre of Yee cell (a_i, b_j, c_k).

In a 2D simulation, J_tot = 0. This means we cannot interpolate in two directions for the Hx and Hz fields, since there is no y-direction to interpolate in. As a result, we only interpolate in the z-direction for Hx and x-direction for Hz when running a two-dimensional simulation.
*/

double MagneticSplitField::interpolate_x_to_centre(int i, int j, int k) {
  if (J_tot==0) {
    // in a 2D simulation, we must interpolate in the z-direction to recover Hx, due to the magnetic-field offsets from the centre
    // determine the interpolation scheme to use
    const InterpolationScheme &scheme = best_scheme(K_tot, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
      interp_data[ind] = xy[k - scheme.number_of_datapoints_to_left + ind][j][i] +
                         xz[k - scheme.number_of_datapoints_to_left + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    return scheme.interpolate(interp_data);
  }
  else {
    // Associations: a = x, b = y, c = z

    // determine the z-direction scheme
    const InterpolationScheme &z_scheme = best_scheme(K_tot, k);
    // determine the y-direction scheme
    const InterpolationScheme &y_scheme = best_scheme(J_tot, j);

    // this data will be passed to the second interpolation scheme
    double data_for_second_scheme[8];
    // this data will hold values for the interpolation in the first interpolation scheme
    double data_for_first_scheme[8];

    // which of z_scheme and y_scheme is WORSE? We will interpolate in this direction SECOND
    if (z_scheme.is_better_than(y_scheme)) {
      // we will be interpolating in the z-direction first, then in y
      for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++) {
        // this is the j-index of the cell we are looking at
        int cell_j = j - y_scheme.number_of_datapoints_to_left + jj;
        // determine the Hx values of cells (i, cell_j, k-z_scheme.index-1) through (i, cell_j, k-z_scheme.index-1+7), and interpolate them via z_scheme
        for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++) {
          // the k-index of the current cell we are looking at (readability, don't need to define this here)
          int cell_k = k - z_scheme.number_of_datapoints_to_left + kk;
          // gather the data for interpolating in the z dimension
          data_for_first_scheme[kk] = xy[cell_k][cell_j][i] + xz[cell_k][cell_j][i];
        }
        // interpolate in z to obtain a value for the Hx field at position (i, cell_j+Dy, k)
        // place this into the appropriate index in the data being passed to the y_scheme
        data_for_second_scheme[jj] = z_scheme.interpolate(data_for_first_scheme);
      }
      // now interpolate in the y-direction to the centre of Yee cell (i,j,k)
      return y_scheme.interpolate(data_for_second_scheme);
    } else {
      // we will be interpolating in the y-direction first, then in z
      for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++) {
        // this is the k-index of the cell we are looking at
        int cell_k = k - z_scheme.number_of_datapoints_to_left + kk;
        // determine the Hx values of cells (i, j - y_scheme.index-1, cell_k) through (i, j - y_scheme.index-1+7, cell_k), and interpolate them via y_scheme
        for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++) {
          // the j-index of the current cell we are looking at (readability, don't need to define this here)
          int cell_j = j - y_scheme.number_of_datapoints_to_left + jj;
          // gather the data for interpolating in the y dimension
          data_for_first_scheme[jj] = xy[cell_k][cell_j][i] + xz[cell_k][cell_j][i];
        }
        // interpolate in y to obtain a value for the Hx field at position (i, j, cell_k+Dz)
        // place this into the appropriate index in the data being passed to the y_scheme
        data_for_second_scheme[kk] = y_scheme.interpolate(data_for_first_scheme);
      }
      // now interpolate in the z-direction to the centre of Yee cell (i,j,k)
      return z_scheme.interpolate(data_for_second_scheme);
    }
  }
}
double MagneticSplitField::interpolate_y_to_centre(int i, int j, int k) {
  // we can always interpolate in two directions for Hy, since even in a 2D simulation the x- and z-directions still exist
  // Associations: a = y, b = z, c = x

  // determine the x-direction scheme
  const InterpolationScheme &x_scheme = best_scheme(I_tot, i);
  // determine the z-direction scheme
  const InterpolationScheme &z_scheme = best_scheme(K_tot, k);

  // this data will be passed to the second interpolation scheme
  double data_for_second_scheme[8];
  // this data will hold values for the interpolation in the first interpolation scheme
  double data_for_first_scheme[8];

  // which of x_scheme and z_scheme is WORSE? We will interpolate in this direction SECOND
  if (z_scheme.is_better_than(x_scheme)) {
    // we will be interpolating in the z-direction first, then in x
    for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++) {
      // this is the i-index of the cell we are looking at
      int cell_i = i - x_scheme.number_of_datapoints_to_left + ii;
      // determine the Hy values of cells (cell_i, j, k-z_scheme.index-1) through (cell_i, j, k-z_scheme.index-1+7), and interpolate them via z_scheme
      for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++) {
        // the k-index of the current cell we are looking at (readability, don't need to define this here)
        int cell_k = k - z_scheme.number_of_datapoints_to_left + kk;
        // gather the data for interpolating in the z dimension
        data_for_first_scheme[kk] = yx[cell_k][j][cell_i] + yz[cell_k][j][cell_i];
      }
      // interpolate in z to obtain a value for the Hy field at position (cell_i+Dx, j, k)
      // place this into the appropriate index in the data being passed to the x_scheme
      data_for_second_scheme[ii] = z_scheme.interpolate(data_for_first_scheme);
    }
    // now interpolate in the x-direction to the centre of Yee cell (i,j,k)
    return x_scheme.interpolate(data_for_second_scheme);
  } else {
    // we will be interpolating in the x-direction first, then in z
    for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++) {
      // this is the k-index of the cell we are looking at
      int cell_k = k - z_scheme.number_of_datapoints_to_left + kk;
      // determine the Hy values of cells (i - x_scheme.index-1, j, cell_k) through (i- x_scheme.index-1+7, j, cell_k), and interpolate them via x_scheme
      for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++) {
        // the i-index of the current cell we are looking at (readability, don't need to define this here)
        int cell_i = i - x_scheme.number_of_datapoints_to_left + ii;
        // gather the data for interpolating in the x dimension
        data_for_first_scheme[ii] = yx[cell_k][j][cell_i] + yz[cell_k][j][cell_i];
      }
      // interpolate in x to obtain a value for the Hy field at position (i, j, cell_k+Dz)
      // place this into the appropriate index in the data being passed to the y_scheme
      data_for_second_scheme[kk] = x_scheme.interpolate(data_for_first_scheme);
    }
    // now interpolate in the z-direction to the centre of Yee cell (i,j,k)
    return z_scheme.interpolate(data_for_second_scheme);
  }
}
double MagneticSplitField::interpolate_z_to_centre(int i, int j, int k) {
  if (J_tot==0) {
    // in a 2D simulation, we must interpolate in the x-direction to recover Hz, due to the magnetic-field offsets from the centre
    // determine the interpolation scheme to use
    const InterpolationScheme &scheme = best_scheme(I_tot, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
      interp_data[ind] = zx[k][j][i - scheme.number_of_datapoints_to_left + ind] +
                         zy[k][j][i - scheme.number_of_datapoints_to_left + ind];
    }

    // now run the interpolation scheme and place the result into the output
    return scheme.interpolate(interp_data);
  }
  else {
    // Associations: a = z, b = x, c = y

    // determine the x-direction scheme
    const InterpolationScheme &x_scheme = best_scheme(I_tot, i);
    // determine the y-direction scheme
    const InterpolationScheme &y_scheme = best_scheme(J_tot, j);

    // this data will be passed to the second interpolation scheme
    double data_for_second_scheme[8];
    // this data will hold values for the interpolation in the first interpolation scheme
    double data_for_first_scheme[8];

    // which of x_scheme and y_scheme is WORSE? We will interpolate in this direction SECOND
    if (y_scheme.is_better_than(x_scheme)) {
      // we will be interpolating in the y-direction first, then in x
      for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++) {
        // this is the i-index of the cell we are looking at
        int cell_i = i - x_scheme.number_of_datapoints_to_left + ii;
        // determine the Hz values of cells (cell_i, j-y_scheme.index-1, k) through (cell_i, j-y_scheme.index-1+7, k), and interpolate them via y_scheme
        for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++) {
          // the j-index of the current cell we are looking at (readability, don't need to define this here)
          int cell_j = j - y_scheme.number_of_datapoints_to_left + jj;
          // gather the data for interpolating in the y dimension
          data_for_first_scheme[jj] = zx[k][cell_j][cell_i] + zy[k][cell_j][cell_i];
        }
        // interpolate in y to obtain a value for the Hz field at position (cell_i+Dx, j, k)
        // place this into the appropriate index in the data being passed to the x_scheme
        data_for_second_scheme[ii] = y_scheme.interpolate(data_for_first_scheme);
      }
      // now interpolate in the x-direction to the centre of Yee cell (i,j,k)
      return x_scheme.interpolate(data_for_second_scheme);
    } else {
      // we will be interpolating in the x-direction first, then in y
      for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++) {
        // this is the j-index of the cell we are looking at
        int cell_j = j - y_scheme.number_of_datapoints_to_left + jj;
        // determine the Hz values of cells (i - x_scheme.index-1, cell_j, k) through (i- x_scheme.index-1+7, cell_j, k), and interpolate them via x_scheme
        for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++) {
          // the i-index of the current cell we are looking at (readability, don't need to define this here)
          int cell_i = i - x_scheme.number_of_datapoints_to_left + ii;
          // gather the data for interpolating in the x dimension
          data_for_first_scheme[ii] = zx[k][cell_j][cell_i] + zy[k][cell_j][cell_i];
        }
        // interpolate in x to obtain a value for the Hz field at position (i, j, cell_k+Dz)
        // place this into the appropriate index in the data being passed to the y_scheme
        data_for_second_scheme[jj] = x_scheme.interpolate(data_for_first_scheme);
      }
      // now interpolate in the y-direction to the centre of Yee cell (i,j,k)
      return y_scheme.interpolate(data_for_second_scheme);
    }
  }
}
