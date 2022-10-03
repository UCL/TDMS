/**
 * @file interpolate_Hfield.h
 * @brief Functions for interpolating the Electric field.
 */
#pragma once
#include <complex.h>

/**
 * @brief Interpolate the x-component of the magnetic field to the centre of Yee cell i,j,k.
 * 
 * Time domain corresponds to interpolation of real-valued fields.
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hx-component at position (i, j+Dy, k+Dz) is associated to the Yee cell (i,j,k).
 * To interpolate Hx to the centre of cell (i,j,k), we must therefore interpolate in the y-direction, then the z-direction, or vice-versa.
 *
 * @param[in] Hxy,Hxz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ,nK Number of Yee cells in the j,k directions respectively
 * @param[out] Hx Interpolated value of the Hx field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int nJ, int nK, double *Hx);
/**
 * @brief Interpolate the x-component of the magnetic field to the centre of Yee cell i,j,k.
 *
 * Frequency domain corresponds to interpolation of complex-valued fields.
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hx-component at position (i, j+Dy, k+Dz) is associated to the Yee cell (i,j,k).
 * To interpolate Hx to the centre of cell (i,j,k), we must therefore interpolate in the y-direction, then the z-direction, or vice-versa.
 *
 * @param[in] Hx Hx component associated to the Yee cells.
 * @param[in] i,j,k Yee cell index
 * @param[in] nJ,nK Number of Yee cells in the j,k directions respectively
 * @param[out] Hx_interp Interpolated value of the Hx field at centre of Yee cell i,j,k
 */
void interpolateFreqDomainHx(std::complex<double> ***Hx, int i, int j, int k, int nJ, int nK, std::complex<double> *Hx_interp);

/**
 * @brief Interpolate the y-component of the magnetic field to the centre of Yee cell i,j,k.
 *
 * Time domain corresponds to interpolation of real-valued fields.
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hy-component at position (i+Dx, j, k+Dz) is associated to the Yee cell (i,j,k).
 * To interpolate Hy to the centre of cell (i,j,k), we must therefore interpolate in the x-direction, then the z-direction, or vice-versa.
 *
 * @param[in] Hyx,Hyz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nI,nK Number of Yee cells in the i,k directions respectively
 * @param[out] Hy Interpolated value of the Hy field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHy(double ***Hxy, double ***Hxz, int i, int j, int k, int nI, int nK, double *Hy);
/**
 * @brief Interpolate the y-component of the magnetic field to the centre of Yee cell i,j,k.
 *
 * Frequency domain corresponds to interpolation of complex-valued fields.
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hy-component at position (i+Dx, j, k+Dz) is associated to the Yee cell (i,j,k).
 * To interpolate Hy to the centre of cell (i,j,k), we must therefore interpolate in the x-direction, then the z-direction, or vice-versa.
 *
 * @param[in] Hy Hy component associated to the Yee cells.
 * @param[in] i,j,k Yee cell index
 * @param[in] nI,nK Number of Yee cells in the i,k directions respectively
 * @param[out] Hy_interp Interpolated value of the Hy field at centre of Yee cell i,j,k
 */
void interpolateFreqDomainHy(std::complex<double> ***Hy, int i, int j, int k, int nI, int nK, std::complex<double> *Hy_interp);

/**
 * @brief Interpolate the z-component of the magnetic field to the centre of Yee cell i,j,k.
 *
 * Time domain corresponds to interpolation of real-valued fields.
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hz-component at position (i+Dx, j+Dy, k) is associated to the Yee cell (i,j,k).
 * To interpolate Hz to the centre of cell (i,j,k), we must therefore interpolate in the x-direction, then the y-direction, or vice-versa.
 *
 * @param[in] Hzx,Hzy Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] nI,nJ Number of Yee cells in the i,k directions respectively
 * @param[out] Hz Interpolated value of the Hz field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHz(double ***Hxy, double ***Hxz, int i, int j, int k, int nI, int nJ, double *Hz);
/**
 * @brief Interpolate the z-component of the magnetic field to the centre of Yee cell i,j,k.
 *
 * Frequency domain corresponds to interpolation of complex-valued fields.
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hz-component at position (i+Dx, j+Dy, k) is associated to the Yee cell (i,j,k).
 * To interpolate Hz to the centre of cell (i,j,k), we must therefore interpolate in the x-direction, then the y-direction, or vice-versa.
 *
 * @param[in] Hz Hz component associated to the Yee cells.
 * @param[in] i,j,k Yee cell index
 * @param[in] nI,nJ Number of Yee cells in the i,j directions respectively
 * @param[out] Hz_interp Interpolated value of the Hz field at centre of Yee cell i,j,k
 */
void interpolateFreqDomainHx(std::complex<double> ***Hx, int i, int j, int k, int nJ, int nK, std::complex<double> *Hz_interp);

/**
 * @brief Interpolate the H-field to the centre of Yee cell i,j,k.
 *
 * Time domain corresponds to interpolation of real-valued fields.
 *
 * @param[in] Hxy,Hxz,Hyx,Hyz,Hzx,Hzy Split components of the Yee cell
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] nI,nJ,nK Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Hx,Hy,Hz Interpolated values of the x,y,z (respectively) field component
 */
void interpolateTimeDomainHField(double ***Hxy, double ***Hxz, double ***Hyx,
                                 double ***Hyz, double ***Hzx, double ***Hzy,
                                 int i, int j, int k, int nI, int nJ, int nK,
                                 double *Hx, double *Hy, double *Hz);
/**
 * @brief Interpolate the H-field to the centre of Yee cell i,j,k.
 *
 * Frequency domain corresponds to interpolation of complex-valued fields.
 *
 * @param[in] Hx,Hy,Hz Field components associated to the Yee cells
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] nI,nJ,nK Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Hx_interp,Hy_interp,Hz_interp Interpolated values of the x,y,z (respectively) field component
 */
void interpolateFreqDomainHField(std::complex<double> ***Hx, std::complex<double> ***Hy, std::complex<double> ***Hz,
                                 int i, int j, int k, int nI, int nJ, int nK,
                                 std::complex<double> *Hx_interp, std::complex<double> *Hy_interp, std::complex<double> *Hz_interp);