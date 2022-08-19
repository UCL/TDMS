#include "interpolate_Hfield.h"
#include "interpolation_methods.h"
#include <stdlib.h>

void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx) {
    // first, we need to determine the non-x dimension that we will interpolate into, in order to bring the fields parallel to the Yee cell centre
    // then, we will interpolate in the x-direction

    // determine the schemes available to us in the y- and z-directions
    const interpScheme &z_scheme = best_interp_scheme(K, k), &y_scheme = best_interp_scheme(J, j);
    
    // determine the scheme that we will be using for the x-direction
    const interpScheme &x_scheme = best_interp_scheme(I, i);

    // now we need to assemble the data to pass to x_scheme
    // like in the E-field case, we might reserve more memory than we needed if we happen to be using cubic interpolation, but this isn't the end of the world
    double pass_to_x_scheme[8];
    // likewise, this array will hold the data for the interpolation in the non-x directions
    double other_dim_data[8];

    // for each cell between i - (x_scheme.index+1) + x_scheme.first_nonzero_coeff to i - (x_scheme.index) + x_scheme.last_nonzero_coeff, we need to interpolate in the 2nd (coordinate) dimension
    // we determine this dimension here to avoid an if statement at every value of the loop variable
    // note that because we are only "moving cells" in the x-direction, whichever of z_scheme and y_scheme is better at one cell will also be the best at every other cell.
    if (z_scheme.is_better_than(y_scheme)) {
        // we will be interpolating in the z-direction first, then in x
        // we need to do this for cells k - (z_scheme.index+1) + z_scheme.first_nonzero_coeff through to k - (z_scheme.index+1) + z_scheme.last_nonzero_coeff
        for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
        {
            for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
            {
                // gather the data for interpolating in the z dimension through cell i - (x_scheme.index + 1) + ii
                other_dim_data[kk] = Hxy[k - (z_scheme.index + 1) + kk][j][i - (x_scheme.index + 1) + ii] + Hxz[k - (z_scheme.index + 1) + kk][j][i - (x_scheme.index + 1) + ii];
            }
            // interpolate in z to obtain a value in cell i - (x_scheme.index+1) + ii for the Hx field
            // place this into the appropriate index in the data being passed to the x_scheme
            pass_to_x_scheme[ii] = z_scheme.interpolate(other_dim_data);
        }
    }
    else {
        // we will be interpolating in the y-direction first, then in x
        // we need to do this for cells j - (y_scheme.index+1) + y_scheme.first_nonzero_coeff through to j - (y_scheme.index+1) + y_scheme.last_nonzero_coeff
        for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
        {
            for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
            {
                // gather the data for interpolating in the y dimension through cell i - (x_scheme.index + 1) + ii
                other_dim_data[jj] = Hxy[k][j - (y_scheme.index + 1) + jj][i - (x_scheme.index + 1) + ii] + Hxz[k][j - (y_scheme.index + 1) + jj][i - (x_scheme.index + 1) + ii];
            }
            // interpolate in y to obtain a value in cell i - (x_scheme.index+1) + ii for the Hx field
            // place this into the appropriate index in the data being passed to the x_scheme
            pass_to_x_scheme[ii] = y_scheme.interpolate(other_dim_data);
        }
    }

    // now interpolate in the x-direction to the centre of the Yell cell i,j,k
    *Hx = x_scheme.interpolate(pass_to_x_scheme);
}

void interpolateTimeDomainHy(double ***Hyx, double ***Hyz, int i, int j, int k, int I, int J, int K, double *Hy)
{
    // first, we need to determine the non-y dimension that we will interpolate into, in order to bring the fields parallel to the Yee cell centre
    // then, we will interpolate in the y-direction

    // determine the schemes available to us in the x- and z-directions
    const interpScheme &z_scheme = best_interp_scheme(K, k), &x_scheme = best_interp_scheme(I, i);

    // determine the scheme that we will be using for the y-direction
    const interpScheme &y_scheme = best_interp_scheme(J, j);

    // now we need to assemble the data to pass to y_scheme
    // like in the E-field case, we might reserve more memory than we needed if we happen to be using cubic interpolation, but this isn't the end of the world
    double pass_to_y_scheme[8];
    // likewise, this array will hold the data for the interpolation in the non-y directions
    double other_dim_data[8];

    // for each cell between j - (y_scheme.index+1) + y_scheme.first_nonzero_coeff to j - (y_scheme.index) + y_scheme.last_nonzero_coeff, we need to interpolate in the 2nd (coordinate) dimension
    // we determine this dimension here to avoid an if statement at every value of the loop variable
    // note that because we are only "moving cells" in the y-direction, whichever of z_scheme and x_scheme is better at one cell will also be the best at every other cell.
    if (z_scheme.is_better_than(x_scheme))
    {
        // we will be interpolating in the z-direction first, then in y
        // we need to do this for cells k - (z_scheme.index+1) + z_scheme.first_nonzero_coeff through to k - (z_scheme.index+1) + z_scheme.last_nonzero_coeff
        for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
        {
            for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
            {
                // gather the data for interpolating in the z dimension through cell j - (y_scheme.index + 1) + jj
                other_dim_data[kk] = Hyx[k - (z_scheme.index + 1) + kk][j - (y_scheme.index + 1) + jj][i] + Hyz[k - (z_scheme.index + 1) + kk][j - (y_scheme.index + 1) + jj][i];
            }
            // interpolate in z to obtain a value in cell j - (y_scheme.index+1) + jj for the Hy field
            // place this into the appropriate index in the data being passed to the y_scheme
            pass_to_y_scheme[jj] = z_scheme.interpolate(other_dim_data);
        }
    }
    else
    {
        // we will be interpolating in the x-direction first, then in y
        // we need to do this for cells i - (x_scheme.index+1) + x_scheme.first_nonzero_coeff through to i - (x_scheme.index+1) + x_scheme.last_nonzero_coeff
        for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
        {
            for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
            {
                // gather the data for interpolating in the y dimension through cell j - (y_scheme.index + 1) + jj
                other_dim_data[ii] = Hyx[k][j - (y_scheme.index + 1) + jj][i - (x_scheme.index + 1) + ii] + Hyz[k][j - (y_scheme.index + 1) + jj][i - (x_scheme.index + 1) + ii];
            }
            // interpolate in y to obtain a value in cell j - (y_scheme.index+1) + jj for the Hy field
            // place this into the appropriate index in the data being passed to the y_scheme
            pass_to_y_scheme[jj] = x_scheme.interpolate(other_dim_data);
        }
    }

    // now interpolate in the x-direction to the centre of the Yell cell i,j,k
    *Hy = y_scheme.interpolate(pass_to_y_scheme);
}

void interpolateTimeDomainHz(double ***Hzx, double ***Hzy, int i, int j, int k, int I, int J, int K, double *Hz)
{
    // first, we need to determine the non-z dimension that we will interpolate into, in order to bring the fields parallel to the Yee cell centre
    // then, we will interpolate in the z-direction

    // determine the schemes available to us in the y- and x-directions
    const interpScheme &x_scheme = best_interp_scheme(I, i), &y_scheme = best_interp_scheme(J, j);

    // determine the scheme that we will be using for the z-direction
    const interpScheme &z_scheme = best_interp_scheme(K, k);

    // now we need to assemble the data to pass to z_scheme
    // like in the E-field case, we might reserve more memory than we needed if we happen to be using cubic interpolation, but this isn't the end of the world
    double pass_to_z_scheme[8];
    // likewise, this array will hold the data for the interpolation in the non-z directions
    double other_dim_data[8];

    // for each cell between k - (z_scheme.index+1) + z_scheme.first_nonzero_coeff to k - (z_scheme.index) + z_scheme.last_nonzero_coeff, we need to interpolate in the 2nd (coordinate) dimension
    // we determine this dimension here to avoid an if statement at every value of the loop variable
    // note that because we are only "moving cells" in the z-direction, whichever of x_scheme and y_scheme is better at one cell will also be the best at every other cell.
    if (x_scheme.is_better_than(y_scheme))
    {
        // we will be interpolating in the x-direction first, then in z
        // we need to do this for cells i - (x_scheme.index+1) + x_scheme.first_nonzero_coeff through to i - (x_scheme.index+1) + x_scheme.last_nonzero_coeff
        for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
        {
            for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
            {
                // gather the data for interpolating in the x dimension through cell k - (z_scheme.index + 1) + kk
                other_dim_data[ii] = Hzx[k - (z_scheme.index + 1) + kk][j][i - (x_scheme.index + 1) + ii] + Hzy[k - (z_scheme.index + 1) + kk][j][i - (x_scheme.index + 1) + ii];
            }
            // interpolate in x to obtain a value in cell k - (z_scheme.index+1) + kk for the Hz field
            // place this into the appropriate index in the data being passed to the z_scheme
            pass_to_z_scheme[kk] = x_scheme.interpolate(other_dim_data);
        }
    }
    else
    {
        // we will be interpolating in the y-direction first, then in z
        // we need to do this for cells j - (y_scheme.index+1) + y_scheme.first_nonzero_coeff through to j - (y_scheme.index+1) + y_scheme.last_nonzero_coeff
        for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
        {
            for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
            {
                // gather the data for interpolating in the y dimension through cell k - (z_scheme.index + 1) + kk
                other_dim_data[jj] = Hzx[k - (z_scheme.index + 1) + kk][j - (y_scheme.index + 1) + jj][i] + Hzy[k - (z_scheme.index + 1) + kk][j - (y_scheme.index + 1) + jj][i];
            }
            // interpolate in y to obtain a value in cell k - (z_scheme.index+1) + kk for the Hz field
            // place this into the appropriate index in the data being passed to the z_scheme
            pass_to_z_scheme[kk] = y_scheme.interpolate(other_dim_data);
        }
    }

    // now interpolate in the x-direction to the centre of the Yell cell i,j,k
    *Hz = z_scheme.interpolate(pass_to_z_scheme);
}

void interpolateTimeDomainHField(double ***Hxy, double ***Hxz, double ***Hyx,
                                 double ***Hyz, double ***Hzx, double ***Hzy,
                                 int i, int j, int k, int I, int J, int K,
                                 double *Hx, double *Hy, double *Hz) {
    interpolateTimeDomainHx(Hxy, Hxz, i, j, k, I, J, K, Hx);
    interpolateTimeDomainHy(Hyx, Hyz, i, j, k, I, J, K, Hy);
    interpolateTimeDomainHz(Hzx, Hzy, i, j, k, I, J, K, Hz);
}