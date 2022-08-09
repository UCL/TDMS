double interp1(double v1, double v2, double v3, double v4);
double interp1(double *v);
double interp2(double v1, double v2, double v3, double v4);
double interp2(double *v);
double interp3(double v1, double v2, double v3, double v4);
double interp3(double *v);
void checkInterpolationPoints(int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int I, int J, int K);

// define aliases for interpolation scheme flags, for readability
enum interp_scheme
{
    BAND_LIMITED_3 = 0, // use bandlimited_interpolation w/ interp position = 3
    BAND_LIMITED_2 = 1,
    BAND_LIMITED_4 = 2,
    BAND_LIMITED_1 = 3,
    BAND_LIMITED_5 = 4,
    BAND_LIMITED_0 = 5,
    BAND_LIMITED_6 = 6,
    CUBIC_INTERP_1 = 7, // use cubic interpolation w/ interp position 1 
    CUBIC_INTERP_2 = 8,
    CUBIC_INTERP_3 = 9
};

interp_scheme determineInterpScheme(int cells_in_direction, int cell_id);

double bandlimited_interpolation(int interp_pos, double *a);