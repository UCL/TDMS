#include "interpolation_methods.h"

#include <stdexcept>
#include <string>

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

void checkInterpolationPoints(int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int nI, int nJ, int nK)
{
    if (i_l < 2)
    {
        throw runtime_error("Interpolation error: i_l too small");
    }
    else if (i_u > nI - 2)
    {
        throw runtime_error("Interpolation error: i_u too large");
    }
    else if (j_l < 2)
    {
        throw runtime_error("Interpolation error: j_l too small");
    }
    else if (j_u > nJ - 2)
    {
        throw runtime_error("Interpolation error: j_u too large");
    }
    else if (k_l < 2)
    {
        throw runtime_error("Interpolation error: k_l too small");
    }
    else if (k_u > nK - 2)
    {
        throw runtime_error("Interpolation error: k_u too large");
    }
}

// THESE WILL BE THE RELEVANT FUNCTIONS ONCE INTERPOLATION SCHEMES ARE COMPLETE AND TESTED

interpScheme::interpScheme(scheme_value val) {

    // set the value field
    priority = val;

    // set values for this method based on the scheme_value passed in
    // these are all hard-coded values - there is no way around a long list of cases!
    switch (val) {
        case BAND_LIMITED_0: 
        {
            scheme_coeffs[0] = 0.389880694606603;
            scheme_coeffs[1] = 0.752494454088346;
            scheme_coeffs[2] = -0.182115867658487;
            scheme_coeffs[3] = 0.046235288061499;
            scheme_coeffs[4] = -0.006777513830539;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 4;
            number_of_datapoints_to_left = 1;
            break; 
        }
        case BAND_LIMITED_1:
        {
            scheme_coeffs[0] = -0.077297572626688;
            scheme_coeffs[1] = 0.570378586429859;
            scheme_coeffs[2] = 0.616613874491358;
            scheme_coeffs[3] = -0.142658093427528;
            scheme_coeffs[4] = 0.039457774230959;
            scheme_coeffs[5] = -0.006777513830539;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 5;
            number_of_datapoints_to_left = 2;
            break;
        }
        case BAND_LIMITED_2:
        {
            scheme_coeffs[0] = 0.025902746569881;
            scheme_coeffs[1] = -0.135880579596988;
            scheme_coeffs[2] = 0.609836360660818;
            scheme_coeffs[3] = 0.609836360660818;
            scheme_coeffs[4] = -0.142658093427528;
            scheme_coeffs[5] = 0.039457774230959;
            scheme_coeffs[6] = -0.006777513830539;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 6;
            number_of_datapoints_to_left = 3;
            break;
        }
        case BAND_LIMITED_3:
        {
            scheme_coeffs[0] = -0.006777513830539;
            scheme_coeffs[1] = 0.039457774230959;
            scheme_coeffs[2] = -0.142658093427528;
            scheme_coeffs[3] = 0.609836360660818;
            scheme_coeffs[4] = 0.609836360660818;
            scheme_coeffs[5] = -0.142658093427528;
            scheme_coeffs[6] = 0.039457774230959;
            scheme_coeffs[7] =  -0.006777513830539;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 7;
            number_of_datapoints_to_left = 4;
            break;
        }
        case BAND_LIMITED_4: 
        {
            scheme_coeffs[1] = -0.006777513830539;
            scheme_coeffs[2] = 0.039457774230959;
            scheme_coeffs[3] = -0.142658093427528;
            scheme_coeffs[4] = 0.609836360660818;
            scheme_coeffs[5] = 0.609836360660818;
            scheme_coeffs[6] = -0.135880579596988;
            scheme_coeffs[7] = 0.025902746569881;
            first_nonzero_coeff = 1;
            last_nonzero_coeff = 7;
            number_of_datapoints_to_left = 5;
            break;
        }
        case BAND_LIMITED_5:
        {
            scheme_coeffs[2] = -0.006777513830539;
            scheme_coeffs[3] = 0.039457774230959;
            scheme_coeffs[4] = -0.142658093427528;
            scheme_coeffs[5] = 0.616613874491358;
            scheme_coeffs[6] = 0.570378586429859;
            scheme_coeffs[7] = -0.077297572626688;
            first_nonzero_coeff = 2;
            last_nonzero_coeff = 7;
            number_of_datapoints_to_left = 6;
            break;
        }
        case BAND_LIMITED_6:
        {
            scheme_coeffs[3] = -0.006777513830539;
            scheme_coeffs[4] = 0.046235288061499;
            scheme_coeffs[5] = -0.182115867658487;
            scheme_coeffs[6] = 0.752494454088346;
            scheme_coeffs[7] = 0.389880694606603;
            first_nonzero_coeff = 3;
            last_nonzero_coeff = 7;
            number_of_datapoints_to_left = 7;
            break;
        }
        case BAND_LIMITED_7:
        {
            scheme_coeffs[3] = 0.006777513830539;
            scheme_coeffs[4] = -0.046235288061499;
            scheme_coeffs[5] = 0.182115867658487;
            scheme_coeffs[6] = -0.752494454088346;
            scheme_coeffs[7] = 1.609553415928240;
            first_nonzero_coeff = 3;
            last_nonzero_coeff = 7;
            number_of_datapoints_to_left = 8;
            break;
        }
        case BAND_LIMITED_CELL_ZERO:
        {
            // by symmetry, we can reflect the setup for BLi to position 7 in order to interpolate to the centre of Yee cell 0 
            scheme_coeffs[0] = 1.609553415928240;
            scheme_coeffs[1] = -0.752494454088346;
            scheme_coeffs[2] = 0.182115867658487;
            scheme_coeffs[3] = -0.046235288061499;
            scheme_coeffs[4] = 0.006777513830539;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 4;
            number_of_datapoints_to_left = 0;
            break;
        }
        case CUBIC_INTERP_FIRST:
        {
            scheme_coeffs[0] = 5. / 16.;
            scheme_coeffs[1] = 15. / 16.;
            scheme_coeffs[2] = -5. / 16.;
            scheme_coeffs[3] = 1. / 16.;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 3;
            number_of_datapoints_to_left = 1;
            break;
        }
        case CUBIC_INTERP_MIDDLE:
        {
            scheme_coeffs[0] = -1. / 16.;
            scheme_coeffs[1] = 9. / 16.;
            scheme_coeffs[2] = 9. / 16.;
            scheme_coeffs[3] = -1. / 16.;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 3;
            number_of_datapoints_to_left = 2;
            break;
        }
        case CUBIC_INTERP_LAST:
        {
            scheme_coeffs[0] = 1. / 16.;
            scheme_coeffs[1] = -5. / 16.;
            scheme_coeffs[2] = 15. / 16.;
            scheme_coeffs[3] = 5. / 16.;
            first_nonzero_coeff = 0;
            last_nonzero_coeff = 3;
            number_of_datapoints_to_left = 3;
            break;
        }
        default:
        {
            // if we cannot identify the scheme from it's value, throw an error
            throw runtime_error("Error: could not assign value " + std::to_string(val) + " to interpolation scheme.\n");
            break;
        }
    }
}

scheme_value interpScheme::get_priority() const {
    return priority;
}

int interpScheme::num_nonzero_coeffs() const {
    return last_nonzero_coeff - first_nonzero_coeff + 1;
}

double interpScheme::interpolate(const double *v, const int offset) const {

    double interp_value = 0.;
    for(int ind=first_nonzero_coeff; ind<=last_nonzero_coeff; ind++) {
        interp_value += scheme_coeffs[ind] * v[ind+offset];
    }
    return interp_value;
}
complex<double> interpScheme::interpolate(const complex<double> *v, const int offset) const {
    complex<double> interp_value = 0.;
    for (int ind = first_nonzero_coeff; ind <= last_nonzero_coeff; ind++) {
        interp_value += scheme_coeffs[ind] * v[ind + offset];
    }
    return interp_value;
}

bool interpScheme::is_better_than(const interpScheme s) const {
    return (priority > s.get_priority());
}

const interpScheme &best_interp_scheme(int cells_in_direction, int cell_id) {

    // interpolation is impossible with fewer than 4 cells in a dimension
    if (cells_in_direction < 4) {
        throw out_of_range("Error: computational domain has fewer than 4 cells in at least 1 dimension, cubic and bandlimited interpolation impossible.\n");
    }
    // Yee cell with index <0 does exist (starts from 0)
    // Current code functionality is to error if we attempt to interpolate to cell 0, because the cubic interpolation found this impossible
    // However, BL_TO_CELL_0 can handle this, provided cells_in_direction >= 8.
    else if (cell_id <= 0) {
        throw out_of_range("Error: Yee cell index <=0 requested (must be >=1).\n");
    }
    // Yee cell with index >= cells_in_direction doesn't exist
    // Again, we can in theory determine the value "here" using BL7, however cubic interpolation cannot do this, and so current code functionality is to throw an error here.
    else if (cell_id >= cells_in_direction) {
        throw out_of_range("Error: Yee cell index beyond maximum number of Yee cells requested.\n");
    }
    else if (cells_in_direction < 8) {
        // we do not have enough cells to use bandlimited interpolation, but can use cubic
        // by definition, cell_id = 1 requires us to use CBFst, cell_id = cells_in_direction-1 CBLast,
        // and everything else CBMid
        if (cell_id == 1) {
            return CBFst;
        }
        else if (cell_id == cells_in_direction-1) {
            return CBLst;
        }
        else {
            return CBMid;
        }
    }
    else {
        // we can apply bandlimited interpolation.
        // unless we are <=3 cells away from either boundary, we will want to use BL3
        if ((cell_id >= 4) && (cell_id <= cells_in_direction - 4)) {
            return BL3;
        }
        else if (cell_id == 1) {
            return BL0;
        }
        else if (cell_id == 2) {
            return BL1;
        }
        else if (cell_id == 3) {
            return BL2;
        }
        else if (cell_id == cells_in_direction - 3) {
            return BL4;
        }
        else if (cell_id == cells_in_direction - 2) {
            return BL5;
        }
        else if (cell_id == cells_in_direction - 1) {
            return BL6;
        }
        else {
            // we somehow got to here, but we should have covered all possible bases above. Return an error
            throw runtime_error("Error: could not identify scheme despite appropriate Yee cell index and number of cells.\n");
        }
    }
}
