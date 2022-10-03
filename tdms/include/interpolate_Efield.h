/**
 * @file interpolate_Efield.h
 * @brief Functions for interpolating the Electric field.
 */
#pragma once

#include <complex.h>

/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell.
 * 
 * Time domain corresponds to interpolating real-valued fields.
 *
 * @param[in] Exy,Exz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nI Number of Yee cells in the x-dimension
 * @param[out] Ex Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEx(double ***Exy, double ***Exz, int i, int j, int k, int nI, double *Ex);
/**
 * @brief Interpolate the Ex field component to the centre of a Yee cell.
 * 
 * Frequency domain corresponds to interpolating complex-valued fields.
 *
 * @param[in] Ex Ex components associated to the Yee cells
 * @param[in] i,j,k Yee cell index
 * @param[in] nI Number of Yee cells in the x-dimension
 * @param[out] Ex_interp Interpolated value of the Ex field at centre of Yee cell i,j,k
 */
void interpolateFreqDomainEx(std::complex<double> ***Ex, int i, int j, int k, int nI, std::complex<double> *Ex_interp);

/**
 * @brief Interpolate the Ey field component to the centre of a Yee cell.
 *
 * Time domain corresponds to interpolating real-valued fields.
 *
 * @param[in] Eyx,Eyz split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ Number of Yee cells in the y-dimension
 * @param[out] Ey Interpolated value of the Ey field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEy(double ***Eyx, double ***Eyz, int i, int j, int k, int nJ, double *Ey);
/**
 * @brief Interpolate the Ey field component to the centre of a Yee cell.
 *
 * Frequency domain corresponds to interpolating complex-valued fields.
 *
 * @param[in] Ey Ey components associated to the Yee cells
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ Number of Yee cells in the y-dimension
 * @param[out] Ey_interp Interpolated value of the Ey field at centre of Yee cell i,j,k
 */
void interpolateFreqDomainEx(std::complex<double> ***Ey, int i, int j, int k, int nJ, std::complex<double> *Ey_interp);

/**
 * @brief Interpolate the Ez field component to the centre of a Yee cell.
 *
 * Time domain corresponds to interpolating real-valued fields.
 *
 * @param[in] Ezx,Ezy split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nK Number of Yee cells in the z-dimension
 * @param[out] Ez Interpolated value of the Ez field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainEz(double ***Ezx, double ***Ezy, int i, int j, int k, int nK, double *Ez);
/**
 * @brief Interpolate the Ez field component to the centre of a Yee cell.
 *
 * Frequency domain corresponds to interpolating complex-valued fields.
 *
 * @param[in] Ez Ez components associated to the Yee cells
 * @param[in] i,j,k Yee cell index
 * @param[in] nK Number of Yee cells in the x-dimension
 * @param[out] Ez_interp Interpolated value of the Ez field at centre of Yee cell i,j,k
 */
void interpolateFreqDomainEx(std::complex<double> ***Ez, int i, int j, int k, int nK, std::complex<double> *Ez_interp);

/**
 * @brief Interpolate the E-field to the centre of Yee cell i,j,k
 * 
 * Time domain corresponds to interpolating real-valued fields.
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
/**
 * @brief Interpolate the E-field to the centre of Yee cell i,j,k
 *
 * Frequency domain corresponds to interpolating complex-valued fields.
 *
 * @param[in] Ex,Ey,Ez Field component values associated to the Yee cells
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] nI,nJ,nK Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Ex_interp,Ey_interp,Ez_interp Interpolated values of the x,y,z (respectively) field component
 */
void interpolateFreqDomainEField(std::complex<double> ***Ex, std::complex<double> ***Ey, std::complex<double> ***Ez,
                                 int i, int j, int k, int nI, int nJ, int nK,
                                 std::complex<double> *Ex_interp, std::complex<double> *Ey_interp, std::complex<double> *Ez_interp);