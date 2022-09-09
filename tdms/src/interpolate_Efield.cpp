#include "interpolate_Efield.h"
#include "interpolation_methods.h"


void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int I, double *Ex) {

    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(I, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // i + scheme.index-1 is the index of the Yee cell that plays the role of v0 in the interpolation
    // note that scheme.index-1 appears here because Yee cells associate field values "to the right" of their centre (IE, the field component is at a position with a higher coordinate value than the Yee cell centre)
    for(int ind=scheme.first_nonzero_coeff; ind<=scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Exy[k][j][i - (scheme.index-1) + ind] + Exz[k][j][i - (scheme.index-1) + ind];
    }

    // now run the interpolation scheme and place the result into the output
    *Ex = scheme.interpolate(interp_data);
}

void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int J, double *Ey)
{
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(J, j);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // j + scheme.index-1 is the index of the Yee cell that plays the role of v0 in the interpolation
    // note that scheme.index-1 appears here because Yee cells associate field values "to the right" of their centre (IE, the field component is at a position with a higher coordinate value than the Yee cell centre)
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++)
    {
        interp_data[ind] = Eyx[k][j - (scheme.index - 1) + ind][i] + Eyz[k][j - (scheme.index - 1) + ind][j];
    }

    // now run the interpolation scheme and place the result into the output
    *Ey = scheme.interpolate(interp_data);
}

void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int K, double *Ez) {

    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(K, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // k + scheme.index-1 is the index of the Yee cell that plays the role of v0 in the interpolation
    // note that scheme.index-1 appears here because Yee cells associate field values "to the right" of their centre (IE, the field component is at a position with a higher coordinate value than the Yee cell centre)
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++)
    {
        interp_data[ind] = Ezx[k - (scheme.index - 1) + ind][j][i] + Ezy[k - (scheme.index - 1) + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ez = scheme.interpolate(interp_data);
}

void interpolateTimeDomainEField(double ***Exy, double ***Exz, double ***Eyx,
                                 double ***Eyz, double ***Ezx, double ***Ezy,
                                 int i, int j, int k, int I, int J, int K,
                                 double *Ex, double *Ey, double *Ez) 
{
    interpolateTimeDomainEx(Exy, Exz, i, j, k, I, Ex);
    interpolateTimeDomainEy(Eyx, Eyz, i, j, k, J, Ey);
    interpolateTimeDomainEz(Ezx, Ezy, i, j, k, K, Ez);
}