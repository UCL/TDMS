# include "interpolation_methods.h"
# include "stdexcept"
# include <string>

using namespace std;

// CUBIC INTERPOLATION FUNCTIONS

/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp1(double v1, double v2, double v3, double v4)
{
    return -1. / 16. * v1 + 9. / 16. * v2 + 9. / 16. * v3 - 1. / 16. * v4;
}
/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp1(double *v) {
    return -1. / 16. * v[0] + 9. / 16. * v[1] + 9. / 16. * v[2] - 1. / 16. * v[3];
}

/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
 */
double interp2(double v1, double v2, double v3, double v4)
{
    return 5. / 16. * v1 + 15. / 16. * v2 - 5. / 16. * v3 + 1. / 16. * v4;
}
/*Use cubic interpolation to interpolate between the last 2 of 4 points
 * v0    v1    v2    v3
 * o     o     o  x  o
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
 * @brief Determines the appropriate interpolation scheme to use
 *
 * @param[in] cells_in_direction The number of Yee cells in the interpolation direction of interest
 * @param[in] cell_id The current ID (in this dimension) of the Yee cell
 * @return interp_scheme {BAND_LIMITED, INTERP1, INTERP2, INTERP3} Indicating the appropriate interpolation scheme
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
    else if (cells_in_direction < 8) {
        // we do not have enough cells to use bandlimited interpolation, but we can use cubic
        // by definition, cell_id = 1 requires us to use interp2, cell_id = 2 interp1, and cell_id = 3 interp3
        switch(cell_id) {
            case 1:
                return CUBIC_INTERP_2;
                break;
            case 2:
                return CUBIC_INTERP_1;
                break;
            case 3:
                return CUBIC_INTERP_3;
                break;
        }
    }
    else {
        // we can apply bandlimited interpolation
        if ((cell_id >= 4) && (cell_id <= cells_in_direction-3)) {return BAND_LIMITED_3;} // best, and most frequent, case
        else if (cell_id==1) {return BAND_LIMITED_0;}
        else if (cell_id == 2) {return BAND_LIMITED_1;}
        else if (cell_id == 3) {return BAND_LIMITED_2;}
        else if (cell_id == cells_in_direction-3) {return BAND_LIMITED_4;}
        else if (cell_id == cells_in_direction-2) {return BAND_LIMITED_5;}
        else {return BAND_LIMITED_6;} //cell_id = cells_in_direction-1
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
        for (int i; i < 8; i++)
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
            // interpolate a[0.5]
            const double b05[5] = {0.389880694606603, 0.752494454088346, -0.182115867658487,
                                   0.046235288061499, -0.006777513830539}; // last 3 terms are 0
            for (int i = 0; i < 5; i++)
            {
                interp_value += a[i+offset] * b05[i];
            }
            break;
        case 1:
            // interpolate a[1.5]
            const double b15[6] = {-0.077297572626688, 0.570378586429859, 0.616613874491358,
                                   -0.142658093427528, 0.039457774230959, -0.006777513830539}; // last 2 terms are 0
            for (int i = 0; i < 6; i++)
            {
                interp_value += a[i+offset] * b15[i];
            }
            break;
        case 2:
            // interpolate a[2.5]
            const double b25[7] = {0.025902746569881, -0.135880579596988, 0.609836360660818,
                                   0.609836360660818, -0.142658093427528, 0.039457774230959,
                                   0.006777513830539};
            for (int i = 0; i < 7; i++)
            {
                interp_value += a[i+offset] * b25[i];
            }
            break;
        case 4:
            // interpolate a[4.5]
            const double b45[7] = {-0.006777513830539, 0.039457774230959, -0.142658093427528,
                                   0.609836360660818, 0.609836360660818, -0.135880579596988,
                                   0.025902746569881}; // first term is zero
            for (int i = 1; i < 8; i++)
            {
                interp_value += a[i+offset] * b45[i - 1];
            }
            break;
        case 5:
            // interpolate a[5.5]
            const double b55[6] = {-0.006777513830539, 0.039457774230959, -0.142658093427528,
                                   0.616613874491358, 0.570378586429859, -0.077297572626688}; // first two terms are zero
            for (int i = 2; i < 8; i++)
            {
                interp_value += a[i+offset] * b55[i - 2];
            }
            break;
        case 6:
            // interpolate a[6.5]
            const double b65[5] = {-0.006777513830539, 0.046235288061499, -0.182115867658487,
                                   0.752494454088346, 0.389880694606603}; // first 3 terms are zero
            for (int i = 3; i < 8; i++)
            {
                interp_value += a[i+offset] * b65[i - 3];
            }
            break;
        case 7:
            // interpolate a[7.5]
            const double b75[5] = {0.006777513830539, -0.046235288061499, 0.182115867658487,
                                   -0.752494454088346, 1.609553415928240}; // first 3 terms are zero
            for (int i = 3; i < 8; i++)
            {
                interp_value += a[i+offset] * b75[i - 3];
            }
            break;
        default:
            throw out_of_range("Requested bandlimited interpolation to position " + to_string((double)interp_pos + 0.5) + ", which is outside limits [0.5,7.5]\n");
            break;
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