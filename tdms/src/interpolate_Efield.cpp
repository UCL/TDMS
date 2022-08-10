#include "interpolate_Efield.h"
#include "interpolation_methods.h"

/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell
 * 
 * @param[in] Exy,Exz split components of the Yee cell 
 * @param[in] i,j,k Yee cell index 
 * @param[in] I Number of Yee cells in the x-dimension
 * @param[out] Ex Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int I, double *Ex) {
    // determine the interpolation scheme to use
    interp_scheme scheme = determineInterpScheme(I, i);
    // prepare input data
    double interp_data_bli[8];
    double interp_data_cubic[4];
    // now fill the interpolation data, depending on whether we are doing cubic or BL interpolation
    // due to how we have setup interp_scheme, we can check whether we are doing BL or cubic interpolation through a comparison
    if (scheme < CUBIC_INTERP_FIRST) {
        // performing BLi, write interpolation data
        for(int ind=0; ind<8; ind++) {
            // we need Yee cell data for cells i - (scheme+1) through i - (scheme+1) + 7, inclusive
            interp_data_bli[i] = Exy[k][j][i - (scheme+1) + ind] + Exz[k][j][i - (scheme+1) + ind];
        }
        // now perform BLi and place the result into the output
        *Ex = bandlimited_interpolation(scheme, interp_data_bli, 0);
    }
    else {
        // performing cubic interpolation, write interpolation data
        for(int ind=0; ind<4; ind++) {
            // we need Yee cell data for cells i - (scheme-7+1) through i - (scheme-7+1) + 3, inclusive
            // scheme -6 is used to save on the extra int computation
            interp_data_cubic[i] = Exy[k][j][i - (scheme - 6) + ind] + Exz[k][j][i - (scheme - 6) + ind];
        }
        // now perform cubic interpolation and place the result into the output
        *Ex = cubic_interpolation(scheme, interp_data_cubic);
    }
}

/**
 * @brief Interpolate the Ey field component to the centre of a Yee cell
 *
 * @param[in] Eyx,Eyz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] J Number of Yee cells in the y-dimension
 * @param[out] Ey Interpolated value of the Ey field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int J, double *Ey)
{
    // determine the interpolation scheme to use
    interp_scheme scheme = determineInterpScheme(J, j);
    // prepare input data
    double interp_data_bli[8];
    double interp_data_cubic[4];
    // now fill the interpolation data, depending on whether we are doing cubic or BL interpolation
    // due to how we have setup interp_scheme, we can check whether we are doing BL or cubic interpolation through a comparison
    if (scheme < CUBIC_INTERP_FIRST)
    {
        // performing BLi, write interpolation data
        for (int ind = 0; ind < 8; ind++)
        {
            // we need Yee cell data for cells j - (scheme+1) through j - (scheme+1) + 7, inclusive
            interp_data_bli[i] = Eyx[k][j - (scheme + 1) + ind][i] + Eyz[k][j - (scheme + 1) + ind][i];
        }
        // now perform BLi and place the result into the output
        *Ey = bandlimited_interpolation(scheme, interp_data_bli, 0);
    }
    else
    {
        // performing cubic interpolation, write interpolation data
        for (int ind = 0; ind < 4; ind++)
        {
            // we need Yee cell data for cells j - (scheme-7+1) through j - (scheme-7+1) + 3, inclusive
            // scheme -6 is used to save on the extra int computation
            interp_data_cubic[i] = Eyx[k][j - (scheme - 6) + ind][i] + Eyz[k][j - (scheme - 6) + ind][i];
        }
        // now perform cubic interpolation and place the result into the output
        *Ey = cubic_interpolation(scheme, interp_data_cubic);
    }
}

/**
 * @brief Interpolate the Ez field component to the centre of a Yee cell
 * 
 * @param[in] Ezx,Ezy split components of the Yee cell 
 * @param[in] i,j,k Yee cell index 
 * @param[in] K Number of Yee cells in the z-dimension
 * @param[out] Ez Interpolated value of the Ez field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int K, double *Ez) {
    // determine the interpolation scheme to use
    interp_scheme scheme = determineInterpScheme(K, k);
    // prepare input data
    double interp_data_bli[8];
    double interp_data_cubic[4];
    // now fill the interpolation data, depending on whether we are doing cubic or BL interpolation
    // due to how we have setup interp_scheme, we can check whether we are doing BL or cubic interpolation through a comparison
    if (scheme < CUBIC_INTERP_FIRST) {
        // performing BLi, write interpolation data
        for(int ind=0; ind<8; ind++) {
            // we need Yee cell data for cells k - (scheme+1) through k - (scheme+1) + 7, inclusive
            interp_data_bli[i] = Ezx[k - (scheme + 1) + ind][j][i] + Ezy[k - (scheme + 1) + ind][j][i];
        }
        // now perform BLi and place the result into the output
        *Ez = bandlimited_interpolation(scheme, interp_data_bli, 0);
    }
    else {
        // performing cubic interpolation, write interpolation data
        for(int ind=0; ind<4; ind++) {
            // we need Yee cell data for cells k - (scheme-7+1) through k - (scheme-7+1) + 3, inclusive
            // scheme -6 is used to save on the extra int computation
            interp_data_cubic[i] = Ezx[k - (scheme - 6) + ind][j][i] + Ezy[k - (scheme - 6) + ind][j][i];
        }
        // now perform cubic interpolation and place the result into the output
        *Ez = cubic_interpolation(scheme, interp_data_cubic);
    }
}

/**
 * @brief Interpolate the E-field to the centre of Yee cell i,j,k
 * 
 * @param[in] Exy,Exz,Eyx,Eyz,Ezx,Ezy Split components of the Yee cell 
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of 
 * @param[in] I,J,K Number of Yee cells in the i,j,k directions (respectively) 
 * @param[out] Ex,Ey,Ez Interpolated values of the x,y,z (respectively) field component
 */
void interpolateTimeDomainEField(double ***Exy, double ***Exz, double ***Eyx,
                                 double ***Eyz, double ***Ezx, double ***Ezy,
                                 int i, int j, int k, int I, int J, int K,
                                 double *Ex, double *Ey, double *Ez) 
{
    interpolateTimeDomainEx(Exy, Exz, i, j, k, I, Ex);
    interpolateTimeDomainEy(Eyx, Eyz, i, j, k, J, Ey);
    interpolateTimeDomainEz(Ezx, Ezy, i, j, k, K, Ez);
}