/**
 * @brief Interpolate the x-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions - first we must interpolate to obtain nearby field values parallel to the centre of the Yee cell, then we must interpolate again in another dimension to the centre of the Yee cell itself.
 *
 * @param[in] Hxy,Hxz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] I,J,K Number of Yee cells in the i,j,k directions respectively
 * @param[out] Hx Interpolated value of the Hx field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx);

/**
 * @brief Interpolate the y-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions - first we must interpolate to obtain nearby field values parallel to the centre of the Yee cell, then we must interpolate again in another dimension to the centre of the Yee cell itself.
 *
 * @param[in] Hyx,Hyz Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] I,J,K Number of Yee cells in the i,j,k directions respectively
 * @param[out] Hy Interpolated value of the Hy field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHy(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx);

/**
 * @brief Interpolate the z-component of the magnetic field to the centre of Yee cell i,j,k
 *
 * Interpolation for the magnetic field must be done across two dimensions - first we must interpolate to obtain nearby field values parallel to the centre of the Yee cell, then we must interpolate again in another dimension to the centre of the Yee cell itself.
 *
 * @param[in] Hzx,Hzy Split components of the Yee cell
 * @param[in] i,j,k Yee cell index
 * @param[in] I,J,K Number of Yee cells in the i,j,k directions respectively
 * @param[out] Hz Interpolated value of the Hz field at centre of Yee cell i,j,k
 */
void interpolateTimeDomainHz(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx);

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
