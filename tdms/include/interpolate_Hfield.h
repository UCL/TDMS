/**
 * @brief Interpolate the x-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hx-component at position (i, j+Dy, k+Dz) is associated to the Yee cell (i,j,k).
 * To interpolate Hx to the centre of cell (i,j,k), we must therefore interpolate in the y-direction, then the z-direction, or vice-versa.
 *
 * @param[in] Hxy,Hxz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] J,K Number of Yee cells in the j,k directions respectively
 * @param[out] Hx Interpolated value of the Hx field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int J, int K, double *Hx);

/**
 * @brief Interpolate the y-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hy-component at position (i+Dx, j, k+Dz) is associated to the Yee cell (i,j,k).
 * To interpolate Hy to the centre of cell (i,j,k), we must therefore interpolate in the x-direction, then the z-direction, or vice-versa.
 *
 * @param[in] Hyx,Hyz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] I,K Number of Yee cells in the i,k directions respectively
 * @param[out] Hy Interpolated value of the Hy field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHy(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int K, double *Hy);

/**
 * @brief Interpolate the z-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions.
 * The stored Hz-component at position (i+Dx, j+Dy, k) is associated to the Yee cell (i,j,k).
 * To interpolate Hz to the centre of cell (i,j,k), we must therefore interpolate in the x-direction, then the y-direction, or vice-versa.
 *
 * @param[in] Hzx,Hzy Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] I,J Number of Yee cells in the i,j directions respectively
 * @param[out] Hz Interpolated value of the Hz field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHz(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, double *Hz);

/**
 * @brief Interpolate the H-field to the centre of Yee cell i,j,k
 *
 * @param[in] Hxy,Hxz,Hyx,Hyz,Hzx,Hzy Split components of the Yee cell
 * @param[in] i,j.k Index of the Yee cell to interpolate to the centre of
 * @param[in] I,J,K Number of Yee cells in the i,j,k directions (respectively)
 * @param[out] Hx,Hy,Hz Interpolated values of the x,y,z (respectively) field component
 */
void interpolateTimeDomainHField(double ***Hxy, double ***Hxz, double ***Hyx,
                                 double ***Hyz, double ***Hzx, double ***Hzy,
                                 int i, int j, int k, int I, int J, int K,
                                 double *Hx, double *Hy, double *Hz);
