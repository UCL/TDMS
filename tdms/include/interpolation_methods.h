// CUBIC INTERPOLATION FUNCTIONS - LEAVE HERE TO PRESERVE interpolate.cpp current functionality prior to testing!

/*Use cubic interpolation to interpolate between the middle 2 of 4 points
 * v0    v1    v2    v3
 * o     o  x   o    o
 */
double interp1(double v1, double v2, double v3, double v4);
/*Use cubic interpolation to interpolate between the middle 2 of 4 points
 * v0    v1    v2    v3
 * o     o  x   o    o
 */
double interp1(double *v);
/*Use cubic interpolation to interpolate between the first 2 of 4 points
 * v0    v1    v2    v3
 * o  x   o     o    o
 */
double interp2(double v1, double v2, double v3, double v4);
/*Use cubic interpolation to interpolate between the first 2 of 4 points
 * v0    v1    v2    v3
 * o  x   o     o    o
 */
double interp2(double *v);
/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp3(double v1, double v2, double v3, double v4);
/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp3(double *v);

/**
 * @brief Checks whether the limits of field extraction are within range of the FDTD grid
 *
 * Since cubic interpolation is being used, it must be ensured that (in any given direction) that the least index used is no smaller than 2, whilst the greatest no larger than the maximum number of cells in that direction - 2.
 *
 * @param i_l Least i index into the FDTD grid to evaluate the field at
 * @param i_u Greatest i index into the FDTD grid to evaluate the field at
 * @param j_l Least j index into the FDTD grid to evaluate the field at
 * @param j_u Greatest j index into the FDTD grid to evaluate the field at
 * @param k_l Least k index into the FDTD grid to evaluate the field at
 * @param k_u Greatest k index into the FDTD grid to evaluate the field at
 * @param I Number of elements in the i direction of the FDTD grid
 * @param J Number of elements in the j direction of the FDTD grid
 * @param K Number of elements in the k direction of the FDTD grid
 *
 * @throws runtime_error In the event that the limits of field extraction are outside the FDTD grid
 */
void checkInterpolationPoints(int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int I, int J, int K);

/**
 * @brief Defines aliases for the flags that will determine the interpolation scheme to be used.
 *
 * The choice of expansion is not arbitrary: for the band-limited cases, the expanded number can be used to tell us the range of Yee indicies we need to extract values from, to then interpolate.
 * Similarly for cubic interpolation: the expanded number, offset by 7, fulfils the same role.
 */
enum interp_scheme
{
    BAND_LIMITED_0 = 0,      // use bandlimited_interpolation w/ interp position = 0
    BAND_LIMITED_1 = 1,      // use bandlimited_interpolation w/ interp position = 1
    BAND_LIMITED_2 = 2,      // use bandlimited_interpolation w/ interp position = 2
    BAND_LIMITED_3 = 3,      // use bandlimited_interpolation w/ interp position = 3 [Preferred method if available]
    BAND_LIMITED_4 = 4,      // use bandlimited_interpolation w/ interp position = 4
    BAND_LIMITED_5 = 5,      // use bandlimited_interpolation w/ interp position = 5
    BAND_LIMITED_6 = 6,      // use bandlimited_interpolation w/ interp position = 6
    CUBIC_INTERP_MIDDLE = 8, // cubic interpolation to middle 2 of 4 points (interp1)
    CUBIC_INTERP_FIRST = 7,  // cubic interpolation to first 2 of 4 points (interp2)
    CUBIC_INTERP_LAST = 9    // cubic interpolation to last 2 of 4 points (interp3)
};

/**
 * @brief Defines aliases for our preference of interpolation scheme to apply.
 *
 * Lower values correspond to higher preference.
 */
enum scheme_value
{
    VAL_BL_3 = 0,  // highest priority is BLi at position 3
    VAL_BL_42 = 1, // BLi at position 4 or 2
    VAL_BL_51 = 2, // BLi at position 5 or 1
    VAL_BL_60 = 3, // BLi at position 6 or 0
    VAL_C_MID = 4, // Cubic at middle position
    VAL_C_OTH = 5  // Cubic at either 1st or last position
};

/**
 * @brief Holds information about the interpolation scheme that is to be applied to a given dataset.
 * 
 */
class interpScheme {
    public:
        // constructors

        /**
         * @brief Default constructor - assumes we can use bandlimited interpolation
         * 
         */
        interpScheme();
        /**
         * @brief Creates an instance of the interpScheme class, by determining the appropriate interpolation scheme to use
         *
         * @param[in] cells_in_direction The number of Yee cells in the interpolation direction of interest
         * @param[in] cell_id The current ID (in this dimension) of the Yee cell
         * @return interp_scheme Indicating the appropriate interpolation scheme to use
         */
        interpScheme(int cells_in_direction, int cell_id);
        /**
         * @brief Implicit conversion to interpScheme by only providing the scheme to use
         * 
         * @param s Interpolation scheme to use
         */
        interpScheme(interp_scheme s);

        // standard fetch methods

        interp_scheme get_scheme();
        scheme_value get_value();

        // comparison methods

        /**
         * @brief Determines if this scheme has higher preference than the scheme s1.
         * 
         * Colloquially, this is the operation "this > s1".
         *
         * @param s1 Interpolation scheme to compare against
         * @return true This is the optimal scheme.
         * @return false s1 is the optimal scheme, or is of equal preference.
         */
        bool is_better_than(interpScheme s1);

        // interpolation method

        /**
         * @brief Interpolate the data v using the interpolation scheme
         * 
         * @param v Data points to interpolate
         * @param offset Read buffer from v[offset] rather than v[0]
         * @return double Interpolated value
         */
        double interpolate(double *v, int offset = 0);

    private:
        // the interpolation scheme to be applied
        interp_scheme scheme;
        // the preference of this interpolation scheme relative to the other schemes
        scheme_value value;

        // private interpolation methods (used internally on call to interpolate)

        /**
         * @brief Performs bandlimited interpolation with 8 sample points, to position i.5
         *
         * Given equidistant sample points a[0],....,a[7], the bandlimited interpolation to the midpoint of a[i] and a[i+1], denoted a[i.5], is
         *  a[i.5] = \sum_{k=0}^7 a[k] * b^{(i.5)}[k],
         * where (for each i), b^{i.5)}[k] is a vector of constant coefficients coded into this function.
         *
         * @param interp_pos The value i in the formula above
         * @param a Equally spaced sample points, must contain (at least 8) elements
         * @param offset {Default 0} a[offset] will be treated as a[0]
         * @return double The interpolated value a[interp_pos.5]
         */
        double bandlimited_interpolation(double *a, int offset = 0);
        /**
         * @brief Performs bandlimited interpolation with 8 sample points, to position i.5
         *
         * Given equidistant sample points a0,....,a7, the bandlimited interpolation to the midpoint of a[i] and a[i+1], denoted a[i.5], is
         *  a[i.5] = \sum_{k=0}^7 a[k] * b^{(i.5)}[k],
         * where (for each i), b^{(i.5)}[k] is a vector of constant coefficients coded into this function.
         *
         * @param interp_pos The value i in the formula above
         * @param a0,a1,a2,a3,a4,a5,a6,a7 Equidistant sample values
         * @return double The interpolated value a[interp_pos.5]
         */
        double bandlimited_interpolation(double a0, double a1, double a2, double a3,
                                         double a4, double a5, double a6, double a7);

        /**
         * @brief Performs cubic interpolation (to the appropriate interpolation point) given 4 data points.
         *
         * Interpolation can be made to 3 different positions between the 4 data points. The scheme determines which position to interpolate to.
         *
         * scheme               v0      v1      v2      v3
         *  CUBIC_FIRST             x
         *  CUBIC_MIDDLE                     x
         *  CUBIC_LAST                              x
         *
         * @param v Equally spaced data points to use for interpolation
         * @param offset [Default 0] Read buffer starting from a[offset] rather than a[0]
         * @return double Interpolated value
         */
        double cubic_interpolation(double *v, int offset = 0);
        /**
         * @brief Performs cubic interpolation (to the appropriate interpolation point) given 4 data points.
         *
         * Interpolation can be made to 3 different positions between the 4 data points. The scheme determines which position to interpolate to.
         *
         * scheme               v0      v1      v2      v3
         *  CUBIC_FIRST             x
         *  CUBIC_MIDDLE                     x
         *  CUBIC_LAST                              x
         *
         * @param v0,v1,v2,v3 Equally spaced data points to use for interpolation
         * @return double Interpolated value
         */
        double cubic_interpolation(double v0, double v1, double v2, double v3);
};