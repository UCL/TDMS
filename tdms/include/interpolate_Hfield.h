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

/**
 * @brief Determines which of two dimensions is optimal to perform interpolation in.
 *
 * Interpolating the magnetic field requires us to first interpolate to obtain points parallel to the centre of the Yee cell, then again in the direction corresponding to the field component we are interested in.
 * This function determines which of the remaining two dimensions, labelled d0 and d1, is optimal for the first set of interpolations.
 *
 * (d0,d1) represents any combination of (x,y), (x,z), or (z,y) dimensions.
 *
 * @param n_cells_d0,n_cells_d1 The maxmium number of Yee cells in dimensions d0 and d1
 * @param cid0,cid1 The component of the Yee cell index in dimensions d0 and d1
 * @return true Dimension d0 should be used.
 * @return false Dimension d1 should be used. In the event of a tie, prefer this dimension.
 */
bool determine_second_dim(int n_cells_d0, int n_cells_d1, int cid0, int cid1);