#include "interpolate_Hfield.h"
#include "interpolation_methods.h"

void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int J, int K, double *Hx)
{
    /* There is a choice in the order of interpolation; either interpolate in y first, then z, or vice-versa.
    For optimal accuracy, we must determine the best interpolation scheme we can use in the y and z directions at cell (i,j,k).
    The WORSE scheme is the scheme that will be applied second, let us assume WLOG that the y-direction has the worse scheme.
    We now need to interpolate Hx in the z-direction, to position (i, cell_j+Dy, k) of cells (i, j-y_scheme.index-1, k) through (i,j-y_scheme.index-1+7, k). In reality, we can also skip over the cells which fall on indices that coincide with 0-coefficents in the y_scheme.
        Let cell_j be a particular index such that
            j-y_scheme.index-1+y_scheme.first_nonzero_coeff <= cell_j <= j-y_scheme.index-1+y_scheme.last_nonzero_coeff.
        We need to run z_scheme to obtain the value of Hx at (i,cell_j+Dy, k), which requires us to pull the values of Hx from cells (i, cell_j, k-z_scheme.index-1) through (i, cell_j, k-z_scheme.index-1+7).
        Again, we can shorten this loop by looking at the non-zero coefficients in the z_scheme.
    Each time we interpolate in z, we place the value into another array.
    We then pass this array to y_scheme to interpolate on.

    The case when y and z have their roles reversed is similar.
    */

    // determine the z-direction scheme
    const interpScheme &z_scheme = best_interp_scheme(K, k);
    // determine the y-direction scheme
    const interpScheme &y_scheme = best_interp_scheme(J, j);

    // this data will be passed to the second interpolation scheme
    double pass_to_second_scheme[8];
    // this data will hold values for the interpolation in the first interpolation scheme
    double data_for_first_scheme[8];

    // which of z_scheme and y_scheme is WORSE? We will interpolate in this direction SECOND
    if (z_scheme.is_better_than(y_scheme))
    {
        // we will be interpolating in the z-direction first, then in y
        // this is the scenario described above
        for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
        {
            // this is the j-index of the cell we are looking at
            int cell_j = j - (y_scheme.index + 1) + jj;
            // determine the Hx values of cells (i, cell_j, k-z_scheme.index-1) through (i, cell_j, k-z_scheme.index-1+7), and interpolate them via z_scheme
            for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
            {
                // the k-index of the current cell we are looking at (readability, don't need to define this here)
                int cell_k = k - (z_scheme.index + 1) + kk;
                // gather the data for interpolating in the z dimension
                data_for_first_scheme[kk] = Hxy[cell_k][cell_j][i] + Hxz[cell_k][cell_j][i];
            }
            // interpolate in z to obtain a value for the Hx field at position (i, cell_j+Dy, k)
            // place this into the appropriate index in the data being passed to the y_scheme
            pass_to_second_scheme[jj] = z_scheme.interpolate(data_for_first_scheme);
        }
        // now interpolate in the y-direction to the centre of Yee cell (i,j,k)
        *Hx = y_scheme.interpolate(pass_to_second_scheme);
    }
    else
    {
        // we will be interpolating in the y-direction first, then in z
        // this is the converse scenario described above
        for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
        {
            // this is the k-index of the cell we are looking at
            int cell_k = k - (z_scheme.index + 1) + kk;
            // determine the Hx values of cells (i, j - y_scheme.index-1, cell_k) through (i, j - y_scheme.index-1+7, cell_k), and interpolate them via y_scheme
            for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
            {
                // the j-index of the current cell we are looking at (readability, don't need to define this here)
                int cell_j = j - (y_scheme.index + 1) + jj;
                // gather the data for interpolating in the y dimension
                data_for_first_scheme[jj] = Hxy[cell_k][cell_j][i] + Hxz[cell_k][cell_j][i];
            }
            // interpolate in y to obtain a value for the Hx field at position (i, j, cell_k+Dz)
            // place this into the appropriate index in the data being passed to the y_scheme
            pass_to_second_scheme[kk] = y_scheme.interpolate(data_for_first_scheme);
        }
        // now interpolate in the z-direction to the centre of Yee cell (i,j,k)
        *Hx = z_scheme.interpolate(pass_to_second_scheme);
    }
}

void interpolateTimeDomainHy(double ***Hyx, double ***Hyz, int i, int j, int k, int I, int K, double *Hy)
{
    /* There is a choice in the order of interpolation; either interpolate in x first, then z, or vice-versa.
    For optimal accuracy, we must determine the best interpolation scheme we can use in the x and z directions at cell (i,j,k).
    The WORSE scheme is the scheme that will be applied second, let us assume WLOG that the x-direction has the worse scheme.
    We now need to interpolate Hy in the z-direction, to position (cell_i+Dx, j, k) of cells (i - x_scheme.index-1, j, k) through (i-x_scheme.index-1+7,j, k). In reality, we can also skip over the cells which fall on indices that coincide with 0-coefficents in the x_scheme.
        Let cell_i be a particular index such that
            i-x_scheme.index-1+x_scheme.first_nonzero_coeff <= cell_i <= i-x_scheme.index-1+x_scheme.last_nonzero_coeff.
        We need to run z_scheme to obtain the value of Hy at (cell_i+Dx,j, k), which requires us to pull the values of Hy from cells (cell_i, j, k-z_scheme.index-1) through (cell_i, j, k-z_scheme.index-1+7).
        Again, we can shorten this loop by looking at the non-zero coefficients in the z_scheme.
    Each time we interpolate in z, we place the value into another array.
    We then pass this array to x_scheme to interpolate on.

    The case when x and z have their roles reversed is similar.
    */
    // determine the x-direction scheme
    const interpScheme &x_scheme = best_interp_scheme(I, i);
    // determine the z-direction scheme
    const interpScheme &z_scheme = best_interp_scheme(K, k);

    // this data will be passed to the second interpolation scheme
    double pass_to_second_scheme[8];
    // this data will hold values for the interpolation in the first interpolation scheme
    double data_for_first_scheme[8];

    // which of x_scheme and z_scheme is WORSE? We will interpolate in this direction SECOND
    if (z_scheme.is_better_than(x_scheme))
    {
        // we will be interpolating in the z-direction first, then in x
        // this is the scenario described above
        for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
        {
            // this is the i-index of the cell we are looking at
            int cell_i = i - (x_scheme.index + 1) + ii;
            // determine the Hy values of cells (cell_i, j, k-z_scheme.index-1) through (cell_i, j, k-z_scheme.index-1+7), and interpolate them via z_scheme
            for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
            {
                // the k-index of the current cell we are looking at (readability, don't need to define this here)
                int cell_k = k - (z_scheme.index + 1) + kk;
                // gather the data for interpolating in the z dimension
                data_for_first_scheme[kk] = Hyx[cell_k][j][cell_i] + Hyz[cell_k][j][cell_i];
            }
            // interpolate in z to obtain a value for the Hy field at position (cell_i+Dx, j, k)
            // place this into the appropriate index in the data being passed to the x_scheme
            pass_to_second_scheme[ii] = z_scheme.interpolate(data_for_first_scheme);
        }
        // now interpolate in the x-direction to the centre of Yee cell (i,j,k)
        *Hy = x_scheme.interpolate(pass_to_second_scheme);
    }
    else
    {
        // we will be interpolating in the x-direction first, then in z
        // this is the converse scenario described above
        for (int kk = z_scheme.first_nonzero_coeff; kk <= z_scheme.last_nonzero_coeff; kk++)
        {
            // this is the k-index of the cell we are looking at
            int cell_k = k - (z_scheme.index + 1) + kk;
            // determine the Hy values of cells (i - x_scheme.index-1, j, cell_k) through (i- x_scheme.index-1+7, j, cell_k), and interpolate them via x_scheme
            for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
            {
                // the i-index of the current cell we are looking at (readability, don't need to define this here)
                int cell_i = i - (x_scheme.index + 1) + ii;
                // gather the data for interpolating in the x dimension
                data_for_first_scheme[ii] = Hyx[cell_k][j][cell_i] + Hyz[cell_k][j][cell_i];
            }
            // interpolate in x to obtain a value for the Hy field at position (i, j, cell_k+Dz)
            // place this into the appropriate index in the data being passed to the y_scheme
            pass_to_second_scheme[kk] = x_scheme.interpolate(data_for_first_scheme);
        }
        // now interpolate in the z-direction to the centre of Yee cell (i,j,k)
        *Hy = z_scheme.interpolate(pass_to_second_scheme);
    }
}

void interpolateTimeDomainHz(double ***Hzx, double ***Hzy, int i, int j, int k, int I, int J, double *Hz)
{
    /* There is a choice in the order of interpolation; either interpolate in x first, then y, or vice-versa.
    For optimal accuracy, we must determine the best interpolation scheme we can use in the x and y directions at cell (i,j,k).
    The WORSE scheme is the scheme that will be applied second, let us assume WLOG that the x-direction has the worse scheme.
    We now need to interpolate Hz in the y-direction, to position (cell_i+Dx, j, k) of cells (i - x_scheme.index-1, j, k) through (i-x_scheme.index-1+7,j, k). In reality, we can also skip over the cells which fall on indices that coincide with 0-coefficents in the x_scheme.
        Let cell_i be a particular index such that
            i-x_scheme.index-1+x_scheme.first_nonzero_coeff <= cell_i <= i-x_scheme.index-1+x_scheme.last_nonzero_coeff.
        We need to run y_scheme to obtain the value of Hz at (cell_i+Dx,j, k), which requires us to pull the values of Hz from cells (cell_i, j-y_scheme.index-1, k) through (cell_i, j-y_scheme.index-1+7, k).
        Again, we can shorten this loop by looking at the non-zero coefficients in the y_scheme.
    Each time we interpolate in y, we place the value into another array.
    We then pass this array to x_scheme to interpolate on.

    The case when x and y have their roles reversed is similar.
    */
    // determine the x-direction scheme
    const interpScheme &x_scheme = best_interp_scheme(I, i);
    // determine the y-direction scheme
    const interpScheme &y_scheme = best_interp_scheme(J, j);

    // this data will be passed to the second interpolation scheme
    double pass_to_second_scheme[8];
    // this data will hold values for the interpolation in the first interpolation scheme
    double data_for_first_scheme[8];

    // which of x_scheme and y_scheme is WORSE? We will interpolate in this direction SECOND
    if (y_scheme.is_better_than(x_scheme))
    {
        // we will be interpolating in the y-direction first, then in x
        // this is the scenario described above
        for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
        {
            // this is the i-index of the cell we are looking at
            int cell_i = i - (x_scheme.index + 1) + ii;
            // determine the Hz values of cells (cell_i, j-y_scheme.index-1, k) through (cell_i, j-y_scheme.index-1+7, k), and interpolate them via y_scheme
            for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
            {
                // the j-index of the current cell we are looking at (readability, don't need to define this here)
                int cell_j = j - (y_scheme.index + 1) + jj;
                // gather the data for interpolating in the y dimension
                data_for_first_scheme[jj] = Hzx[k][cell_j][cell_i] + Hzy[k][cell_j][cell_i];
            }
            // interpolate in y to obtain a value for the Hz field at position (cell_i+Dx, j, k)
            // place this into the appropriate index in the data being passed to the x_scheme
            pass_to_second_scheme[ii] = y_scheme.interpolate(data_for_first_scheme);
        }
        // now interpolate in the x-direction to the centre of Yee cell (i,j,k)
        *Hz = x_scheme.interpolate(pass_to_second_scheme);
    }
    else
    {
        // we will be interpolating in the x-direction first, then in y
        // this is the converse scenario described above
        for (int jj = y_scheme.first_nonzero_coeff; jj <= y_scheme.last_nonzero_coeff; jj++)
        {
            // this is the j-index of the cell we are looking at
            int cell_j = j - (y_scheme.index + 1) + jj;
            // determine the Hz values of cells (i - x_scheme.index-1, cell_j, k) through (i- x_scheme.index-1+7, cell_j, k), and interpolate them via x_scheme
            for (int ii = x_scheme.first_nonzero_coeff; ii <= x_scheme.last_nonzero_coeff; ii++)
            {
                // the i-index of the current cell we are looking at (readability, don't need to define this here)
                int cell_i = i - (x_scheme.index + 1) + ii;
                // gather the data for interpolating in the x dimension
                data_for_first_scheme[ii] = Hzx[k][cell_j][cell_i] + Hzy[k][cell_j][cell_i];
            }
            // interpolate in x to obtain a value for the Hz field at position (i, j, cell_k+Dz)
            // place this into the appropriate index in the data being passed to the y_scheme
            pass_to_second_scheme[jj] = x_scheme.interpolate(data_for_first_scheme);
        }
        // now interpolate in the y-direction to the centre of Yee cell (i,j,k)
        *Hz = y_scheme.interpolate(pass_to_second_scheme);
    }
}

void interpolateTimeDomainHField(double ***Hxy, double ***Hxz, double ***Hyx,
                                 double ***Hyz, double ***Hzx, double ***Hzy,
                                 int i, int j, int k, int I, int J, int K,
                                 double *Hx, double *Hy, double *Hz)
{
    interpolateTimeDomainHx(Hxy, Hxz, i, j, k, J, K, Hx);
    interpolateTimeDomainHy(Hyx, Hyz, i, j, k, I, K, Hy);
    interpolateTimeDomainHz(Hzx, Hzy, i, j, k, I, J, Hz);
}