#include "interpolation_methods.h"

#include <stdexcept>
#include <string>

using namespace std;

InterpolationScheme::InterpolationScheme(scheme_value val) {

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

scheme_value InterpolationScheme::get_priority() const {
    return priority;
}

int InterpolationScheme::num_nonzero_coeffs() const {
    return last_nonzero_coeff - first_nonzero_coeff + 1;
}

bool InterpolationScheme::is_better_than(const InterpolationScheme s) const {
    return (priority > s.get_priority());
}

const InterpolationScheme &best_scheme(int datapts_in_direction, int interpolation_position) {

    // interpolation is impossible with fewer than 4 datapoints in a dimension
    if (datapts_in_direction < 4) {
        throw out_of_range("Error: domain axis has <4 datapoints, cubic and bandlimited interpolation impossible.\n");
    }
    else if (datapts_in_direction < 8) {
        // we are restricted to cubic interpolation
        if (interpolation_position <= 0 || interpolation_position >= datapts_in_direction) {
            throw out_of_range("Error: Cubic interpolation impossible to position " + to_string(interpolation_position) + "\n");
        }
        else {
            // determine cubic interpolation scheme to use
            if (interpolation_position == 1) {
                return CBFst;
            }
            else if (interpolation_position == datapts_in_direction - 1) {
                return CBLst;
            }
            else {
                return CBMid;
            }
        }
    }
    else if (interpolation_position < 0 || interpolation_position > datapts_in_direction) {
        // cannot interpolate to here using BLi, throw error
        throw out_of_range("Error: BLi interpolation impossible to position " + to_string(interpolation_position) + "\n");
    }
    else {
        // safe to use BLi, figure out which scheme we need
        if (interpolation_position == 0) {
            return BL_TO_CELL_0;
        }
        else if (interpolation_position == datapts_in_direction) {
            return BL7;
        }
        else if (interpolation_position == 1) {
            return BL0;
        }
        else if (interpolation_position == 2) {
            return BL1;
        }
        else if (interpolation_position == 3) {
            return BL2;
        }
        else if (interpolation_position == datapts_in_direction - 3) {
            return BL4;
        }
        else if (interpolation_position == datapts_in_direction - 2) {
            return BL5;
        }
        else if (interpolation_position == datapts_in_direction - 1) {
            return BL6;
        }
        else {
            // we have 4 datapoints either side of where we want to interpolate to
            return BL3;
        }
    }
    // if we get to here we have, somehow, not returned a scheme. Raise an error
    throw runtime_error("Error: could not identify scheme for unknown reasons, diagnose further.\n");
}
