/**
 * @file interpolate_Efield.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Template functions for interpolating the electric field.
 */
#pragma once

#include <complex>

#include "interpolation_methods.h"

/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell, passing in split field components.
 *
 * @param[in] Exy,Exz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nI Number of Yee cells in the x-dimension
 * @param[out] Ex Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
template <typename T>
void interpolateSplitFieldEx(T ***Exy, T ***Exz, int i, int j, int k, int nI, T *Ex) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nI, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Exy[k][j][i - scheme.number_of_datapoints_to_left + ind] + Exz[k][j][i - scheme.number_of_datapoints_to_left + ind];
    }

    // now run the interpolation scheme and place the result into the output
    *Ex = scheme.interpolate(interp_data);
}
/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell.
 *
 * @param[in] Ex Ex components associated to the Yee cells
 * @param[in] i,j,k Yee cell index
 * @param[in] nI Number of Yee cells in the x-dimension
 * @param[out] Ex_interp Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
template <typename T>
void interpolateEx(T ***Ex, int i, int j, int k, int nI, T *Ex_interp) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nI, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ex[k][j][i - scheme.number_of_datapoints_to_left + ind];
    }
    // now run the interpolation scheme and place the result into the output
    *Ex_interp = scheme.interpolate(interp_data);
}

/**
 * @brief Interpolate the Ey field component to the centre of a Yee cell, passing in split field components.
 *
 * @param[in] Eyx,Eyz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ Number of Yee cells in the y-dimension
 * @param[out] Ey Interpolated value of the Ey field at centre of Yee cell i,j,k
 */
template <typename T>
void interpolateSplitFieldEy(T ***Eyx, T ***Eyz, int i, int j, int k, int nJ, T *Ey) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nJ, j);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Eyx[k][j - scheme.number_of_datapoints_to_left + ind][i] + Eyz[k][j - scheme.number_of_datapoints_to_left + ind][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ey = scheme.interpolate(interp_data);
}
/**
 * @brief Interpolate the Ey field component to the centre of a Yee cell.
 *
 * @param[in] Ey Ey components associated to the Yee cells
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ Number of Yee cells in the y-dimension
 * @param[out] Ey_interp Interpolated value of the Ey field at centre of Yee cell i,j,k
 */
template <typename T>
void interpolateEy(T ***Ey, int i, int j, int k, int nJ, T *Ey_interp) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nJ, j);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ey[k][j - scheme.number_of_datapoints_to_left + ind][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ey_interp = scheme.interpolate(interp_data);
}

/**
 * @brief Interpolate the Ez field component to the centre of a Yee cell, passing in split field components.
 *
 * @param[in] Ezx,Ezy split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nK Number of Yee cells in the z-dimension
 * @param[out] Ez Interpolated value of the Ez field at centre of Yee cell i,j,k
 */
template <typename T>
void interpolateSplitFieldEz(T ***Ezx, T ***Ezy, int i, int j, int k, int nK, T *Ez) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nK, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ezx[k - scheme.number_of_datapoints_to_left + ind][j][i] + Ezy[k - scheme.number_of_datapoints_to_left + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ez = scheme.interpolate(interp_data);
}
/**
 * @brief Interpolate the Ez field component to the centre of a Yee cell.
 *
 * @param[in] Ez Ez components associated to the Yee cells
 * @param[in] i,j,k Yee cell index
 * @param[in] nK Number of Yee cells in the x-dimension
 * @param[out] Ez_interp Interpolated value of the Ez field at centre of Yee cell i,j,k
 */
template <typename T>
void interpolateEz(T ***Ez, int i, int j, int k, int nK, T *Ez_interp) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nK, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ez[k - scheme.number_of_datapoints_to_left + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ez_interp = scheme.interpolate(interp_data);
}

/**
 * @brief Interpolate the E-field to the centre of Yee cell (i,j,k); passing in split field components.
 *
 * @param[in] Exy,Exz,Eyx,Eyz,Ezx,Ezy Split components of the Yee cell
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] nI,nJ,nK Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Ex,Ey,Ez Interpolated values of the x,y,z (respectively) field component
 */
template <typename Tx, typename Ty, typename Tz>
void interpolateSplitFieldE(Tx ***Exy, Tx ***Exz,
                            Ty ***Eyx, Ty ***Eyz,
                            Tz ***Ezx, Tz ***Ezy,
                            int i, int j, int k, int nI, int nJ, int nK,
                            Tx *Ex, Ty *Ey, Tz *Ez) {
    interpolateSplitFieldEx(Exy, Exz, i, j, k, nI, Ex);
    interpolateSplitFieldEy(Eyx, Eyz, i, j, k, nJ, Ey);
    interpolateSplitFieldEz(Ezx, Ezy, i, j, k, nK, Ez);
}
/**
 * @brief Interpolate the E-field to the centre of Yee cell i,j,k
 *
 * @param[in] Ex,Ey,Ez Field component values associated to the Yee cells
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] nI,nJ,nK Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Ex_interp,Ey_interp,Ez_interp Interpolated values of the x,y,z (respectively) field component
 */
template <typename Tx, typename Ty, typename Tz>
void interpolateE(Tx ***Ex, Ty ***Ey, Tz ***Ez,
                  int i, int j, int k, int nI, int nJ, int nK,
                  Tx *Ex_interp, Ty *Ey_interp, Tz *Ez_interp) {
    interpolateEx(Ex, i, j, k, nI, Ex_interp);
    interpolateEy(Ey, i, j, k, nJ, Ey_interp);
    interpolateEz(Ez, i, j, k, nK, Ez_interp);
}