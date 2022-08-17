# include "interpolation_methods.h"
# include "stdexcept"
# include <string>

using namespace std;

// CUBIC INTERPOLATION FUNCTIONS

/*Use cubic interpolation to interpolate between the middle 2 of 4 points
 * v0    v1    v2    v3
 * o     o  x   o    o
 */
double interp1(double v1, double v2, double v3, double v4)
{
    return -1. / 16. * v1 + 9. / 16. * v2 + 9. / 16. * v3 - 1. / 16. * v4;
}
/*Use cubic interpolation to interpolate between the middle 2 of 4 points
 * v0    v1    v2    v3
 * o     o  x   o    o
 */
double interp1(double *v) {
    return -1. / 16. * v[0] + 9. / 16. * v[1] + 9. / 16. * v[2] - 1. / 16. * v[3];
}
/*Use cubic interpolation to interpolate between the first 2 of 4 points
 * v0    v1    v2    v3
 * o  x   o     o    o
 */
double interp2(double v1, double v2, double v3, double v4)
{
    return 5. / 16. * v1 + 15. / 16. * v2 - 5. / 16. * v3 + 1. / 16. * v4;
}
/*Use cubic interpolation to interpolate between the first 2 of 4 points
 * v0    v1    v2    v3
 * o  x   o     o    o
 */
double interp2(double *v)
{
    return 5. / 16. * v[0] + 15. / 16. * v[1] - 5. / 16. * v[2] + 1. / 16. * v[3];
}
/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp3(double v1, double v2, double v3, double v4)
{
    return 1. / 16. * v1 - 5. / 16. * v2 + 15. / 16. * v3 + 5. / 16. * v4;
}
/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp3(double *v)
{
    return 1. / 16. * v[0] - 5. / 16. * v[1] + 15. / 16. * v[2] + 5. / 16. * v[3];
}

/**
 * @brief Performs cubic interpolation to interp_pos, given 4 data points.
 * 
 * Use cubic interpolation to interpolate to the middle of a given pair of points.
 * The input data in the array v, and the position of the interpolated value depending on interp_pos, can be visualised as follows:
 * 
 * interp_pos       v0      v1      v2      v3
 *  0                   x
 *  1                           x
 *  2                                   x
 * 
 * @param interp_pos The position to interpolate to
 * @param v Equally spaced data points to use for interpolation
 * @return double Interpolated value
 */
double cubic_interpolation(int interp_pos, double *v) {
    // determine which interpolation scheme is to be used based off interp_pos
    switch (interp_pos)
    {
    case 0:
        return interp2(v);
        break;
    case 1:
        return interp1(v);
        break;
    case 2:
        return interp3(v);
        break;
    default:
        throw out_of_range("Cubic interpolation cannot take place before 1st point, or after last\n");
        break;
    }
}

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
void checkInterpolationPoints(int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int I, int J, int K)
{
    if (i_l < 2)
    {
        throw runtime_error("Interpolation error: i_l too small");
    }
    else if (i_u > I - 2)
    {
        throw runtime_error("Interpolation error: i_u too large");
    }
    else if (j_l < 2)
    {
        throw runtime_error("Interpolation error: j_l too small");
    }
    else if (j_u > J - 2)
    {
        throw runtime_error("Interpolation error: j_u too large");
    }
    else if (k_l < 2)
    {
        throw runtime_error("Interpolation error: k_l too small");
    }
    else if (k_u > K - 2)
    {
        throw runtime_error("Interpolation error: k_u too large");
    }
}

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
double bandlimited_interpolation(int interp_pos, double *a, int offset)
{
    // in the vast majority of cases, we expect to be well outside the PML and thus in this case
    // so it is more efficient for us to check via an if statement before matching to a switch
    double interp_value = 0.;
    if (interp_pos == 3)
    {
        const double b45[8] = {-0.006777513830539, 0.039457774230959, -0.142658093427528, 0.609836360660818,
                               0.609836360660818, -0.142658093427528, 0.039457774230959, -0.006777513830539};
        for (int i = 0; i < 8; i++)
        {
            interp_value += a[i+offset] * b45[i];
        }
    }
    else
    {
        // we want to interpolate a point that is not at position 3.5
        switch (interp_pos)
        {
        case 0:
            {
            // interpolate a[0.5]
            const double b05[5] = {0.389880694606603, 0.752494454088346, -0.182115867658487,
                                   0.046235288061499, -0.006777513830539}; // last 3 terms are 0
            for (int i = 0; i < 5; i++)
            {
                interp_value += a[i+offset] * b05[i];
            }
            break;
            }
        case 1:
            {
            // interpolate a[1.5]
            const double b15[6] = {-0.077297572626688, 0.570378586429859, 0.616613874491358,
                                   -0.142658093427528, 0.039457774230959, -0.006777513830539}; // last 2 terms are 0
            for (int i = 0; i < 6; i++)
            {
                interp_value += a[i+offset] * b15[i];
            }
            break;
            }
        case 2:
            {
            // interpolate a[2.5]
            const double b25[7] = {0.025902746569881, -0.135880579596988, 0.609836360660818,
                                   0.609836360660818, -0.142658093427528, 0.039457774230959,
                                   0.006777513830539};
            for (int i = 0; i < 7; i++)
            {
                interp_value += a[i+offset] * b25[i];
            }
            break;
            }
        case 4:
            {
            // interpolate a[4.5]
            const double b45[7] = {-0.006777513830539, 0.039457774230959, -0.142658093427528,
                                   0.609836360660818, 0.609836360660818, -0.135880579596988,
                                   0.025902746569881}; // first term is zero
            for (int i = 1; i < 8; i++)
            {
                interp_value += a[i+offset] * b45[i - 1];
            }
            break;
            }
        case 5:
            {
            // interpolate a[5.5]
            const double b55[6] = {-0.006777513830539, 0.039457774230959, -0.142658093427528,
                                   0.616613874491358, 0.570378586429859, -0.077297572626688}; // first two terms are zero
            for (int i = 2; i < 8; i++)
            {
                interp_value += a[i+offset] * b55[i - 2];
            }
            break;
            }
        case 6:
            {
            // interpolate a[6.5]
            const double b65[5] = {-0.006777513830539, 0.046235288061499, -0.182115867658487,
                                   0.752494454088346, 0.389880694606603}; // first 3 terms are zero
            for (int i = 3; i < 8; i++)
            {
                interp_value += a[i+offset] * b65[i - 3];
            }
            break;
            }
        case 7:
            {
            // interpolate a[7.5]
            const double b75[5] = {0.006777513830539, -0.046235288061499, 0.182115867658487,
                                   -0.752494454088346, 1.609553415928240}; // first 3 terms are zero
            for (int i = 3; i < 8; i++)
            {
                interp_value += a[i+offset] * b75[i - 3];
            }
            break;
            }
        default:
            {
            throw out_of_range("Requested bandlimited interpolation to position " + to_string((double)interp_pos + 0.5) + ", which is outside limits [0.5,7.5]\n");
            break;
            }
        } // end switch(interp_pos)
    }     // else
    return interp_value;
}

/**
 * @brief Performs bandlimited interpolation with 8 sample points, to position i.5
 *
 * Given equidistant sample points a0,....,a[7, the bandlimited interpolation to the midpoint of a[i] and a[i+1], denoted a[i.5], is
 *  a[i.5] = \sum_{k=0}^7 a[k] * b^{(i.5)}[k],
 * where (for each i), b^{i.5)}[k] is a vector of constant coefficients coded into this function.
 *
 * @param interp_pos The value i in the formula above
 * @param a0,a1,a2,a3,a4,a5,a6,a7 Equidistant sample values
 * @return double The interpolated value a[interp_pos.5]
 */
double bandlimited_interpolation(int interp_pos, double a0, double a1, double a2, double a3,
                                 double a4, double a5, double a6, double a7) 
{
    double a[8] = {a0, a1, a2, a3, a4, a5, a6, a7};
    return bandlimited_interpolation(interp_pos, a, 0);
}

/**
 * @brief Determines the appropriate interpolation scheme to use
 *
 * @param[in] cells_in_direction The number of Yee cells in the interpolation direction of interest
 * @param[in] cell_id The current ID (in this dimension) of the Yee cell
 * @return interp_scheme Indicating the appropriate interpolation scheme to use
 */
interp_scheme determineInterpScheme(int cells_in_direction, int cell_id)
{
    // interpolation is impossible if the total number of cells in this direction is <4
    if (cells_in_direction < 4)
    {
        throw out_of_range("Error: computational domain has fewer than 4 cells in at least 1 dimension, cubic and band-limited interpolation impossible.\n");
    }
    // Yee cell with index <0 doesn't exist (indexing starts from 1), and Yee cell 0 cannot undergo interpolation to the centre
    else if (cell_id <= 0)
    {
        throw out_of_range("Error: Yee cell index <=0 requested (must be >=1).\n");
    }
    // Yee cell with index >=cells_in_direction doesn't exist
    else if (cell_id >= cells_in_direction)
    {
        throw out_of_range("Error: requested Yee cell index beyond maximum number of Yee cells.\n");
    }
    else if (cells_in_direction < 8)
    {
        // we do not have enough cells to use bandlimited interpolation, but we can use cubic
        // by definition, cell_id = 1 requires us to use interp2, cell_id = cells_in_direction-1 interp3,
        // and everything else interp2
        if (cell_id == 1)
        {
            return CUBIC_INTERP_FIRST;
        }
        else if (cell_id == cells_in_direction - 1)
        {
            return CUBIC_INTERP_LAST;
        }
        else
        {
            return CUBIC_INTERP_MIDDLE;
        }
    }
    else
    {
        // we can apply bandlimited interpolation
        if ((cell_id >= 4) && (cell_id <= cells_in_direction - 3))
        {
            return BAND_LIMITED_3;
        } // best, and most frequent, case
        else if (cell_id == 1)
        {
            return BAND_LIMITED_0;
        }
        else if (cell_id == 2)
        {
            return BAND_LIMITED_1;
        }
        else if (cell_id == 3)
        {
            return BAND_LIMITED_2;
        }
        else if (cell_id == cells_in_direction - 3)
        {
            return BAND_LIMITED_4;
        }
        else if (cell_id == cells_in_direction - 2)
        {
            return BAND_LIMITED_5;
        }
        else
        {
            return BAND_LIMITED_6;
        } // cell_id = cells_in_direction-1
    }
}

/**
 * @brief Given two interpolation schemes s0 and s1, performs evaluation of the expression "s0 > s1", in terms of optimality of the scheme.
 *
 * @param s0,s1 interp_schemes
 * @return true If s0 is the optimal scheme.
 * @return false If s1 is the optimal scheme. Also returned in the event of a tie.
 */
bool better_scheme(interp_scheme s0, interp_scheme s1) {
    /* Optimality is as follows:
    BLi @ pos 3.5
    BLi @ pos 2.5 or 4.5
    BLi @ pos 1.5 or 5.5
    BLi @ pos 0.5 or 6.5
    Cubic @ middle
    Cubic @ first or last
    */
   int s0_rank, s1_rank;
   switch (s0)
   {
    case BAND_LIMITED_0: {s0_rank = 4; break;}
    case BAND_LIMITED_1: {s0_rank = 3; break;}
    case BAND_LIMITED_2: {s0_rank = 2; break;}
    case BAND_LIMITED_3: {s0_rank = 1; break;}
    case BAND_LIMITED_4: {s0_rank = 2; break;}
    case BAND_LIMITED_5: {s0_rank = 3; break;}
    case BAND_LIMITED_6: {s0_rank = 4; break;}
    case CUBIC_INTERP_MIDDLE: {s0_rank = 5; break;}
    case CUBIC_INTERP_FIRST: {s0_rank = 6; break;}
    case CUBIC_INTERP_LAST: {s0_rank = 6; break;}
   }
   switch (s1)
   {
    case BAND_LIMITED_0: {s1_rank = 4; break;}
    case BAND_LIMITED_1: {s1_rank = 3; break;}
    case BAND_LIMITED_2: {s1_rank = 2; break;}
    case BAND_LIMITED_3: {s1_rank = 1; break;}
    case BAND_LIMITED_4: {s1_rank = 2; break;}
    case BAND_LIMITED_5: {s1_rank = 3; break;}
    case BAND_LIMITED_6: {s1_rank = 4; break;}
    case CUBIC_INTERP_MIDDLE: {s1_rank = 5; break;}
    case CUBIC_INTERP_FIRST: {s1_rank = 6; break;}
    case CUBIC_INTERP_LAST: {s1_rank = 6; break;}
   }
   return (s0_rank <= s1_rank);
}