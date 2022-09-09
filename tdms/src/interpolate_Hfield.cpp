#include "interpolate_Hfield.h"
#include "interpolation_methods.h"
#include <stdlib.h>

using namespace std;

/**
 * @brief Determines which of two dimensions is optimal to perform interpolation in.
 *
 * Interpolating the magnetic field requires us to first interpolate to obtain points parallel to the centre of the Yee cell, then again in the direction corresponding to the field component we are interested in.
 * This function determines which of the remaining two dimensions, labelled d0 and d1, is optimal for the first set of interpolations.
 *
 * (d0,d1) represents any combination of (x,y), (x,z), or (z,y) dimensions.
 *
 * @param n_cells_d0,n_cells_d1 The maximum number of Yee cells in dimensions d0 and d1
 * @param cid0,cid1 The component of the Yee cell index in dimensions d0 and d1
 * @return true Dimension d0 should be used.
 * @return false Dimension d1 should be used. In the event of a tie, prefer this dimension.
 */
bool determine_second_dim(int n_cells_d0, int n_cells_d1, int cid0, int cid1) {
    // if we are doing a 2D simulation, one of the dimensions might be zero
    // in which case, return the other
    if (n_cells_d0 == 0) return false;
    else if (n_cells_d1 == 0) return true;
    else {
        // so we have a 3D simulation
        interp_scheme d0_scheme = determineInterpScheme(n_cells_d0, cid0);
        interp_scheme d1_scheme = determineInterpScheme(n_cells_d1, cid1);
        // we now need to determine which scheme is the better scheme
        return better_scheme(d0_scheme, d1_scheme);
    }
}

/**
 * @brief Interpolate the x-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions - first we must interpolate to obtain nearby field values parallel to the centre of the Yee cell, then we must interpolate again in another dimension to the centre of the Yee cell itself.
 *
 * @param[in] Hxy,Hxz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] I,J,K Number of Yee cells in the i,j,k directions respectively
 * @param[out] Hx Interpolated value of the Hx field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx) {
    // which interpolation scheme do we want to use in the x-direction?
    interp_scheme x_scheme = determineInterpScheme(I, i);

    // this interpolation scheme will determine how many times we need to interpolate in the other dimension, and the range of x-Yee cells to consider
    int i_low, i_high, d_low, d_high;
    // we will also need to store the result of the interpolation in the first dimension
    double *first_interp_vals;
    double *other_dim_cell_data;

    if (better_scheme(x_scheme, CUBIC_INTERP_MIDDLE)) {
        // if using a BLi scheme in x,
        // we need interpolated data for cells i - (scheme+1) through i - (scheme+1) + 7, inclusive
        i_low = i - (x_scheme+1); i_high = i - (x_scheme+1)+7;
        // assign memory that we need to store these values
        first_interp_vals = (double *)malloc(8*sizeof(double));
    }
    else {
        // we're using a cubic scheme in x
        // we need interpolated data for cells i - (scheme-7+1) through i - (scheme-7+1) + 3, inclusive
        i_low = i - (x_scheme-6); i_high = i - (x_scheme-6) + 3;
        // assign memory that we need to store these values
        first_interp_vals = (double *)malloc(4*sizeof(double));
    }

    // now we need to determine which of the y- and z-directions is optimal to perform the first interpolation in
    // note: preference for avoiding the y-dimension because of possible 2D simulations
    if (determine_second_dim(J, K, j, k)) {
        // prefer y-dimension over z-dimension
        interp_scheme y_scheme = determineInterpScheme(J, j);
        if (better_scheme(y_scheme, CUBIC_INTERP_MIDDLE))
        {
            // using a BLi scheme in y, assign memory we need
            other_dim_cell_data = (double *)malloc(8*sizeof(double));

            // Each time we interpolate in y, we need data for cells j - (scheme+1) through j - (scheme+1) + 7, inclusive
            d_low = j - (y_scheme + 1);
            d_high = j - (y_scheme + 1) + 7;
            for (int ii=i_low; ii<=i_high; ii++) {
                // populate other_dim_cell_data with the data we need from the y-cells
                for(int jj=d_low; jj<=d_high; jj++) {
                    other_dim_cell_data[jj-d_low] = Hxy[k][jj][ii] + Hxz[k][jj][ii];
                }
                first_interp_vals[ii-i_low] = bandlimited_interpolation(y_scheme, other_dim_cell_data);
            }
        }
        else
        {
            // we're using a cubic scheme in y, assign memory we need
            other_dim_cell_data = (double *)malloc(4*sizeof(double));

            // Each time we interpolate in y, we need interpolated data for cells j - (scheme-7+1) through j - (scheme-7+1) + 3, inclusive
            d_low = j - (y_scheme - 6);
            d_high = j - (y_scheme - 6) + 3;
            for (int ii = i_low; ii <= i_high; ii++)
            {
                // populate other_dim_cell_data with the data we need from the y-cells
                for (int jj = d_low; jj <= d_high; jj++)
                {
                    other_dim_cell_data[jj - d_low] = Hxy[k][jj][ii] + Hxz[k][jj][ii];
                }
                first_interp_vals[ii - i_low] = cubic_interpolation(y_scheme-7, other_dim_cell_data);
            }
        }
    }
    else {
        // prefer z-dimension over y-dimension
        interp_scheme z_scheme = determineInterpScheme(K, k);
        if (better_scheme(z_scheme, CUBIC_INTERP_MIDDLE))
        {
            // using a BLi scheme in z, assign memory we need
            other_dim_cell_data = (double *)malloc(8 * sizeof(double));

            // Each time we interpolate in y, we need data for cells j - (scheme+1) through j - (scheme+1) + 7, inclusive
            d_low = j - (z_scheme + 1);
            d_high = j - (z_scheme + 1) + 7;
            for (int ii = i_low; ii <= i_high; ii++)
            {
                // populate other_dim_cell_data with the data we need from the y-cells
                for (int kk = d_low; kk <= d_high; kk++)
                {
                    other_dim_cell_data[kk - d_low] = Hxy[kk][j][ii] + Hxz[kk][j][ii];
                }
                first_interp_vals[ii - i_low] = bandlimited_interpolation(z_scheme, other_dim_cell_data);
            }
        }
        else
        {
            // we're using a cubic scheme in y, assign memory we need
            other_dim_cell_data = (double *)malloc(4 * sizeof(double));

            // Each time we interpolate in y, we need interpolated data for cells j - (scheme-7+1) through j - (scheme-7+1) + 3, inclusive
            d_low = j - (z_scheme - 6);
            d_high = j - (z_scheme - 6) + 3;
            for (int ii = i_low; ii <= i_high; ii++)
            {
                // populate other_dim_cell_data with the data we need from the y-cells
                for (int kk = d_low; kk <= d_high; kk++)
                {
                    other_dim_cell_data[kk - d_low] = Hxy[kk][j][ii] + Hxz[kk][j][ii];
                }
                first_interp_vals[ii - i_low] = cubic_interpolation(z_scheme - 7, other_dim_cell_data);
            }
        }
    }

    // now interpolate to the centre of the Yee cell we are about.
    // NOTE: we could do without this if statement by nesting the previous code block in the first
    // if (better_scheme(x_scheme, CUBIC_INTERP_MIDDLE)) - statement of this kind, but that makes for a hideous looking function, and this is already disgusting as it is
    if (better_scheme(x_scheme, CUBIC_INTERP_MIDDLE)) {
        // using BLi in x
        *Hx = bandlimited_interpolation(x_scheme, first_interp_vals);
    }
    else {
        // using cubic interpolation in x
        *Hx = cubic_interpolation(x_scheme-7, first_interp_vals);
    }
}