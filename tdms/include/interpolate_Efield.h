/**
 * @file interpolate_Efield.h
 * @brief Functions for interpolating the Electric field.
 */
#pragma once

/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell
 *
 * @param[in] Exy,Exz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nI Number of Yee cells in the x-dimension
 * @param[out] Ex Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int nI, double *Ex);

/**
 * @brief Interpolate the Ey field component to the centre of a Yee cell
 *
 * @param[in] Eyx,Eyz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ Number of Yee cells in the y-dimension
 * @param[out] Ey Interpolated value of the Ey field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int nJ, double *Ey);

/**
 * @brief Interpolate the Ez field component to the centre of a Yee cell
 *
 * @param[in] Ezx,Ezy split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nK Number of Yee cells in the z-dimension
 * @param[out] Ez Interpolated value of the Ez field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int nK, double *Ez);

/**
 * @brief Interpolate the E-field to the centre of Yee cell i,j,k
 *
 * @param[in] Exy,Exz,Eyx,Eyz,Ezx,Ezy Split components of the Yee cell
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] nI,nJ,nK Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Ex,Ey,Ez Interpolated values of the x,y,z (respectively) field component
 */
void interpolateTimeDomainEField(double ***Exy, double ***Exz, double ***Eyx,
                                 double ***Eyz, double ***Ezx, double ***Ezy, 
                                 int i, int j, int k, int nI, int nJ, int nK, 
                                 double *Ex, double *Ey, double *Ez);