// cubic interpolation functions
double interp1(double v1, double v2, double v3, double v4);
double interp1(double *v);
double interp2(double v1, double v2, double v3, double v4);
double interp2(double *v);
double interp3(double v1, double v2, double v3, double v4);
double interp3(double *v);

double cubic_interpolation(int interp_pos, double *v);

// error checking (whether interpolation can be performed)

void checkInterpolationPoints(int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int I, int J, int K);

// bandlimited interpolation functions

double bandlimited_interpolation(int interp_pos, double *a, int offset = 0);
double bandlimited_interpolation(int interp_pos, double a0, double a1, double a2, double a3,
                                 double a4, double a5, double a6, double a7);

/**
 * @brief Defines aliases for the flags that will determine the interpolation scheme to be used.
 * 
 * The choice of expansion is not arbitrary: for the band-limited cases, the expanded number can be used to tell us the range of Yee indicies we need to extract values from, to then interpolate.
 * Similarly for cubic interpolation: the expanded number, offset by 7, fulfils the same role.
 */
enum interp_scheme
{
    BAND_LIMITED_0 = 0, // use bandlimited_interpolation w/ interp position = 0
    BAND_LIMITED_1 = 1, // use bandlimited_interpolation w/ interp position = 1
    BAND_LIMITED_2 = 2, // use bandlimited_interpolation w/ interp position = 2
    BAND_LIMITED_3 = 3, // use bandlimited_interpolation w/ interp position = 3 [Preferred method if available]
    BAND_LIMITED_4 = 4, // use bandlimited_interpolation w/ interp position = 4
    BAND_LIMITED_5 = 5, // use bandlimited_interpolation w/ interp position = 5
    BAND_LIMITED_6 = 6, // use bandlimited_interpolation w/ interp position = 6
    CUBIC_INTERP_MIDDLE = 8, // cubic interpolation to middle 2 of 4 points (interp1)
    CUBIC_INTERP_FIRST = 7, // cubic interpolation to first 2 of 4 points (interp2)
    CUBIC_INTERP_LAST = 9 // cubic interpolation to last 2 of 4 points (interp3)
};

interp_scheme determineInterpScheme(int cells_in_direction, int cell_id);
bool better_scheme(interp_scheme s0, interp_scheme s1);