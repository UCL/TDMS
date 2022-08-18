#include "interpolate_Efield.h"
#include "interpolation_methods.h"


void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int I, double *Ex) {

    // determine the interpolation scheme to use
    interpScheme scheme_to_use(I, i);
    interp_scheme scheme = scheme_to_use.get_scheme();
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data, depending on whether we are doing cubic or BL interpolation
    // due to how we have setup interp_scheme, we can check whether we are doing BL or cubic interpolation through a comparison
    if (scheme_to_use.is_better_than(CUBIC_INTERP_MIDDLE)) {
        // performing BLi, write interpolation data
        for(int ind=0; ind<8; ind++) {
            // we need Yee cell data for cells i - (scheme+1) through i - (scheme+1) + 7, inclusive
            interp_data[i] = Exy[k][j][i - (scheme+1) + ind] + Exz[k][j][i - (scheme+1) + ind];
        }
    }
    else {
        // performing cubic interpolation, write interpolation data
        for(int ind=0; ind<4; ind++) {
            // we need Yee cell data for cells i - (scheme-7+1) through i - (scheme-7+1) + 3, inclusive
            // scheme -6 is used to save on the extra int computation
            interp_data[i] = Exy[k][j][i - (scheme - 6) + ind] + Exz[k][j][i - (scheme - 6) + ind];
        }
    }
    // now run the interpolation scheme and place the result into the output
    *Ex = scheme_to_use.interpolate(interp_data);
}

void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int J, double *Ey)
{

    // determine the interpolation scheme to use
    interpScheme scheme_to_use(J, j);
    interp_scheme scheme = scheme_to_use.get_scheme();
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data, depending on whether we are doing cubic or BL interpolation
    // due to how we have setup interp_scheme, we can check whether we are doing BL or cubic interpolation through a comparison
    if (scheme_to_use.is_better_than(CUBIC_INTERP_MIDDLE)) {
        // performing BLi, write interpolation data
        for (int ind = 0; ind < 8; ind++)
        {
            // we need Yee cell data for cells j - (scheme+1) through j - (scheme+1) + 7, inclusive
            interp_data[i] = Eyx[k][j - (scheme + 1) + ind][i] + Eyz[k][j - (scheme + 1) + ind][i];
        }
    }
    else
    {
        // performing cubic interpolation, write interpolation data
        for (int ind = 0; ind < 4; ind++)
        {
            // we need Yee cell data for cells j - (scheme-7+1) through j - (scheme-7+1) + 3, inclusive
            // scheme -6 is used to save on the extra int computation
            interp_data[i] = Eyx[k][j - (scheme - 6) + ind][i] + Eyz[k][j - (scheme - 6) + ind][i];
        }
    }
    // now run the interpolation scheme and place the result into the output
    *Ey = scheme_to_use.interpolate(interp_data);
}

void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int K, double *Ez) {

    // determine the interpolation scheme to use
    interpScheme scheme_to_use(K, k);
    interp_scheme scheme = scheme_to_use.get_scheme();
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data, depending on whether we are doing cubic or BL interpolation
    // due to how we have setup interp_scheme, we can check whether we are doing BL or cubic interpolation through a comparison
    if (scheme_to_use.is_better_than(CUBIC_INTERP_MIDDLE)) {
        // performing BLi, write interpolation data
        for(int ind=0; ind<8; ind++) {
            // we need Yee cell data for cells k - (scheme+1) through k - (scheme+1) + 7, inclusive
            interp_data[i] = Ezx[k - (scheme + 1) + ind][j][i] + Ezy[k - (scheme + 1) + ind][j][i];
        }
    }
    else {
        // performing cubic interpolation, write interpolation data
        for(int ind=0; ind<4; ind++) {
            // we need Yee cell data for cells k - (scheme-7+1) through k - (scheme-7+1) + 3, inclusive
            // scheme -6 is used to save on the extra int computation
            interp_data[i] = Ezx[k - (scheme - 6) + ind][j][i] + Ezy[k - (scheme - 6) + ind][j][i];
        }
    }
    // now run the interpolation scheme and place the result into the output
    *Ez = scheme_to_use.interpolate(interp_data);
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