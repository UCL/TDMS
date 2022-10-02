#include "interpolate_Efield.h"
#include "interpolation_methods.h"


void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int nI, double *Ex) {

    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nI, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for(int ind=scheme.first_nonzero_coeff; ind<=scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Exy[k][j][i - scheme.number_of_datapoints_to_left + ind] + Exz[k][j][i - scheme.number_of_datapoints_to_left + ind];
    }

    // now run the interpolation scheme and place the result into the output
    *Ex = scheme.interpolate(interp_data);
}

void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int nJ, double *Ey)
{
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nJ, j);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++)
    {
        interp_data[ind] = Eyx[k][j - scheme.number_of_datapoints_to_left + ind][i] + Eyz[k][j - scheme.number_of_datapoints_to_left + ind][j];
    }

    // now run the interpolation scheme and place the result into the output
    *Ey = scheme.interpolate(interp_data);
}

void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int nK, double *Ez) {

    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nK, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++)
    {
        interp_data[ind] = Ezx[k - scheme.number_of_datapoints_to_left + ind][j][i] + Ezy[k - scheme.number_of_datapoints_to_left + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ez = scheme.interpolate(interp_data);
}

void interpolateTimeDomainEField(double ***Exy, double ***Exz, double ***Eyx,
                                 double ***Eyz, double ***Ezx, double ***Ezy,
                                 int i, int j, int k, int nI, int nJ, int nK,
                                 double *Ex, double *Ey, double *Ez) 
{
    interpolateTimeDomainEx(Exy, Exz, i, j, k, nI, Ex);
    interpolateTimeDomainEy(Eyx, Eyz, i, j, k, nJ, Ey);
    interpolateTimeDomainEz(Ezx, Ezy, i, j, k, nK, Ez);
}
