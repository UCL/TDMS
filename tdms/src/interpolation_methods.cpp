# include "interpolation_methods.h"
# include "stdexcept"
# include <string>

using namespace std;

// CUBIC INTERPOLATION FUNCTIONS - LEAVE HERE TO PRESERVE interpolate.cpp current functionality prior to testing!

double interp1(double v1, double v2, double v3, double v4)
{
    return -1. / 16. * v1 + 9. / 16. * v2 + 9. / 16. * v3 - 1. / 16. * v4;
}
double interp1(double *v) {
    return -1. / 16. * v[0] + 9. / 16. * v[1] + 9. / 16. * v[2] - 1. / 16. * v[3];
}
double interp2(double v1, double v2, double v3, double v4)
{
    return 5. / 16. * v1 + 15. / 16. * v2 - 5. / 16. * v3 + 1. / 16. * v4;
}
double interp2(double *v)
{
    return 5. / 16. * v[0] + 15. / 16. * v[1] - 5. / 16. * v[2] + 1. / 16. * v[3];
}
double interp3(double v1, double v2, double v3, double v4)
{
    return 1. / 16. * v1 - 5. / 16. * v2 + 15. / 16. * v3 + 5. / 16. * v4;
}
double interp3(double *v)
{
    return 1. / 16. * v[0] - 5. / 16. * v[1] + 15. / 16. * v[2] + 5. / 16. * v[3];
}

void checkInterpolationPoints(int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int I, int J, int K)
{
    if (i_l < 2)
    {
        throw out_of_range("Interpolation error: i_l too small");
    }
    else if (i_u > I - 2)
    {
        throw out_of_range("Interpolation error: i_u too large");
    }
    else if (j_l < 2)
    {
        throw out_of_range("Interpolation error: j_l too small");
    }
    else if (j_u > J - 2)
    {
        throw out_of_range("Interpolation error: j_u too large");
    }
    else if (k_l < 2)
    {
        throw out_of_range("Interpolation error: k_l too small");
    }
    else if (k_u > K - 2)
    {
        throw out_of_range("Interpolation error: k_u too large");
    }
    else if (i_l > i_u) 
    {
        throw runtime_error("Interpolation error: lower cell index greater than upper (i)");
    }
    else if (j_l > j_u && J!=0) // in a 2D simulation (J==0) we allow this beahviour
    {
        throw runtime_error("Interpolation error: lower cell index greater than upper (j)");
    }
    else if (k_l > k_u) 
    {
        throw runtime_error("Interpolation error: lower cell index greater than upper (k)");
    }
}

// THESE WILL BE THE RELEVANT FUNCTIONS ONCE INTERPOLATION SCHEMES ARE COMPLETE AND TESTED

interpScheme::interpScheme() {
    scheme = BAND_LIMITED_3;
    value = VAL_BL_3;
}
interpScheme::interpScheme(int cells_in_direction, int cell_id)
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
            scheme = CUBIC_INTERP_FIRST;
            value = VAL_C_OTH;
        }
        else if (cell_id == (cells_in_direction - 1))
        {
            scheme = CUBIC_INTERP_LAST;
            value = VAL_C_OTH;
        }
        else
        {
            scheme = CUBIC_INTERP_MIDDLE;
            value = VAL_C_MID;
        }
    }
    else
    {
        // we can apply bandlimited interpolation
        if ((cell_id >= 4) && (cell_id <= cells_in_direction - 4))
        {
            scheme = BAND_LIMITED_3;
            value = VAL_BL_3;
        } // best, and most frequent, case
        else if (cell_id == 1)
        {
            scheme = BAND_LIMITED_0;
            value = VAL_BL_60;
        }
        else if (cell_id == 2)
        {
            scheme = BAND_LIMITED_1;
            value = VAL_BL_51;
        }
        else if (cell_id == 3)
        {
            scheme = BAND_LIMITED_2;
            value = VAL_BL_42;
        }
        else if (cell_id == cells_in_direction - 3)
        {
            scheme = BAND_LIMITED_4;
            value = VAL_BL_42;
        }
        else if (cell_id == cells_in_direction - 2)
        {
            scheme = BAND_LIMITED_5;
            value = VAL_BL_51;
        }
        else // cell_id = cells_in_direction-1
        {
            scheme = BAND_LIMITED_6;
            value = VAL_BL_60;
        } 
    }
}
interpScheme::interpScheme(interp_scheme s) {
    scheme = s;
    switch (scheme) {
        case BAND_LIMITED_0:
            value = VAL_BL_60;
            break;
        case BAND_LIMITED_1:
            value = VAL_BL_51;
            break;
        case BAND_LIMITED_2:
            value = VAL_BL_42;
            break;
        case BAND_LIMITED_3:
            value = VAL_BL_3;
            break;
        case BAND_LIMITED_4:
            value = VAL_BL_42;
            break;
        case BAND_LIMITED_5:
            value = VAL_BL_51;
            break;
        case BAND_LIMITED_6:
            value = VAL_BL_60;
            break;
        case CUBIC_INTERP_MIDDLE:
            value = VAL_C_MID;
            break;
        default: // is either cubic-first or cubic-last
            value = VAL_C_OTH;
            break;
    }
}

interp_scheme interpScheme::get_scheme() {
    return scheme;
}

scheme_value interpScheme::get_value() {
    return value;
}

bool interpScheme::is_better_than(interpScheme s1) {
    // lower value == higher preference!
    return ( value < s1.get_value() );
}

double interpScheme::interpolate(double *v, int offset) {
    if (is_better_than(CUBIC_INTERP_MIDDLE)) {
        return bandlimited_interpolation(v, offset);
    }
    else {
        return cubic_interpolation(v, offset);
    }
}

double interpScheme::bandlimited_interpolation(double *a, int offset)
{
    // in the vast majority of cases, we expect to be well outside the PML and thus in this case
    // so it is more efficient for us to check via an if statement before matching to a switch
    double interp_value = 0.;
    if (scheme == BAND_LIMITED_3)
    {
        const double b45[8] = {-0.006777513830539, 0.039457774230959, -0.142658093427528, 0.609836360660818,
                               0.609836360660818, -0.142658093427528, 0.039457774230959, -0.006777513830539};
        for (int i = 0; i < 8; i++)
        {
            interp_value += a[i + offset] * b45[i];
        }
    }
    else
    {
        // we want to interpolate a point that is not at position 3.5
        switch (scheme)
        {
        case 0:
        {
            // interpolate a[0.5]
            const double b05[5] = {0.389880694606603, 0.752494454088346, -0.182115867658487,
                                   0.046235288061499, -0.006777513830539}; // last 3 terms are 0
            for (int i = 0; i < 5; i++)
            {
                interp_value += a[i + offset] * b05[i];
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
                interp_value += a[i + offset] * b15[i];
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
                interp_value += a[i + offset] * b25[i];
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
                interp_value += a[i + offset] * b45[i - 1];
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
                interp_value += a[i + offset] * b55[i - 2];
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
                interp_value += a[i + offset] * b65[i - 3];
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
                interp_value += a[i + offset] * b75[i - 3];
            }
            break;
        }
        default:
        {
            throw out_of_range("Requested bandlimited interpolation to position " + to_string((double)scheme + 0.5) + ", which is outside limits [0.5,7.5]\n");
            break;
        }
        } // end switch(interp_pos)
    }     // else
    return interp_value;
}
double interpScheme::bandlimited_interpolation(double a0, double a1, double a2, double a3,
                                               double a4, double a5, double a6, double a7)
{
    double a[8] = {a0, a1, a2, a3, a4, a5, a6, a7};
    return bandlimited_interpolation(a, 0);
}

double interpScheme::cubic_interpolation(double *v, int offset)
{
    return cubic_interpolation(v[offset], v[offset+1], v[offset+2], v[offset+3]);
}
double interpScheme::cubic_interpolation(double v0, double v1, double v2, double v3) {
    // determine which interpolation scheme is to be used based off scheme
    switch (scheme)
    {
    case CUBIC_INTERP_FIRST:
        return 5. / 16. * v0 + 15. / 16. * v1 - 5. / 16. * v2 + 1. / 16. * v3;
        break;
    case CUBIC_INTERP_MIDDLE:
        return -1. / 16. * v0 + 9. / 16. * v1 + 9. / 16. * v2 - 1. / 16. * v3;
        break;
    case CUBIC_INTERP_LAST:
        return 1. / 16. * v0 - 5. / 16. * v1 + 15. / 16. * v2 + 5. / 16. * v3;
        break;
    default: // should never reach here unless we've called this on a BLi scheme
        throw out_of_range("Cubic interpolation cannot take place before 1st point, or after last\n");
        break;
    }
}
