/**
 * @file interpolation_methods.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief InterpScheme class methods and supporting functions
 *
 * Non InterpScheme methods are required to preserve functionality whilst testing new schemes
 */
#pragma once

#include <complex>

#include "globals.h"

/**
 * @brief Defines our order of preference for the use of the various schemes.
 *
 * There should never be an instance in which we wish to use BLi to position 7 - this will take us to a point OUTSIDE the computational domain. However for completion purposes (and if we find a use for it), it is included.
 *
 * MODIFICATIONS TO THE ALIASED INTS WILL CHANGE THE ORDER OF SCHEME PREFERENCE!
 *
 * For band-limited schemes, we have 8 equidistant datapoints and can interpolate to the midpoint of any pair of consecutive points, or "half a point spacing" to the right of the final data point.
 * These "interpolation positions" (interp positions) are marked with "o" below:
 * v0   v1   v2   v3   v4   v5   v6   v7
 *    o    o    o    o    o    o    o    o
 * We label these interp positions with the index of the point to the left: a data point at interp position 0 is effectively a data point "v0.5", being at the midpoint of v0 and v1, for example.
 */
enum scheme_value
{
    BAND_LIMITED_0 = 4,             // use bandlimited interpolation w/ interp position = 0
    BAND_LIMITED_1 = 6,             // use bandlimited interpolation w/ interp position = 1
    BAND_LIMITED_2 = 8,             // use bandlimited interpolation w/ interp position = 2
    BAND_LIMITED_3 = 9,             // use bandlimited interpolation w/ interp position = 3 [Preferred method if available]
    BAND_LIMITED_4 = 7,             // use bandlimited interpolation w/ interp position = 4
    BAND_LIMITED_5 = 5,             // use bandlimited interpolation w/ interp position = 5
    BAND_LIMITED_6 = 3,             // use bandlimited interpolation w/ interp position = 6
    BAND_LIMITED_7 = 2,             // use bandlimited interpolation w/ interp position = 7
    BAND_LIMITED_CELL_ZERO = 1,     // use bandlimited interpolation to interpolate to the centre of Yee cell 0
    CUBIC_INTERP_MIDDLE = 0,        // cubic interpolation to middle 2 of 4 points (interp1)
    CUBIC_INTERP_FIRST = -1,        // cubic interpolation to first 2 of 4 points (interp2)
    CUBIC_INTERP_LAST = -2          // cubic interpolation to last 2 of 4 points (interp3)
};

class InterpolationScheme {
    private:
        // the "preference" or "value" of applying this scheme. It may be better to apply another scheme with a higher priority.
        scheme_value priority;

        // the constants that will be used in the interpolation scheme.
        double scheme_coeffs[8];
    public:
        /**
         * @brief Construct a new interp Scheme object, by providing the scheme value
         *
         * @param value A value associtated to one of the possible schemes
         */
        InterpolationScheme(scheme_value value);

        /* FETCH METHODS */

        /**
         * @brief Get the value object
         *
         * @return scheme_value
         */
        scheme_value get_priority() const;

        /* END FETCH METHODS */

        /* The number of datapoints "to the left" of where we are planning to interpolate to.
        For example, interpolating to interpolation position 0 requires 1 datapoint at a position before the location of the interpolation, whilst interpolation position 6 requires there to be 7.
        */
        int number_of_datapoints_to_left;

        // cubic and BLi schemes use different numbers of coefficients. To avoid switches, we store these variables.
        int first_nonzero_coeff, last_nonzero_coeff;

        /**
         * @brief Compute the number of non-zero coefficients in the interpolation scheme
         *
         * @return int Number of non-zero coefficients in the interpolation scheme
         */
        int num_nonzero_coeffs() const;

        /**
         * @brief Executes the interpolation scheme on the data provided
         *
         * The interpolation schemes are all of the form
         * interpolated_value = \sum_{i=0}^{7} scheme_coeffs[i] * v[i],
         * so provided that the coefficients have been set correctly in construction (and the data gathered appropriately), we can run the same for loop for each interpolation scheme.
         *
         * For slight speedup, the actual sum performed loops over those i such that
         * 0 <= first_nonzero_coeff <= i <= last_nonzero_coeff <= 7.
         *
         * @param v Sample datapoints to use in interpolation; v[0] should be the first of 8 values
         * @param offset [Default 0] Read buffer from v[offset] rather than v[0]
         * @return double Interpolated value
         */
        template<typename T>
        T interpolate(const T *v, const int offset = 0) const {
          T interp_value = 0.;
          for (int ind = first_nonzero_coeff; ind <= last_nonzero_coeff; ind++) {
            interp_value += scheme_coeffs[ind] * v[ind + offset];
          }
          return interp_value;
        };

        /**
         * @brief Determines whether another interpScheme has greater value than this one
         *
         * @param s The other interpScheme to compare against
         * @return true This scheme has greater value
         * @return false This scheme has lesser, or equal, value to s
         */
        bool is_better_than(const InterpolationScheme s) const;
};

/* Constant members of the interpScheme class */
const InterpolationScheme BL0 = InterpolationScheme(BAND_LIMITED_0);
const InterpolationScheme BL1 = InterpolationScheme(BAND_LIMITED_1);
const InterpolationScheme BL2 = InterpolationScheme(BAND_LIMITED_2);
const InterpolationScheme BL3 = InterpolationScheme(BAND_LIMITED_3);
const InterpolationScheme BL4 = InterpolationScheme(BAND_LIMITED_4);
const InterpolationScheme BL5 = InterpolationScheme(BAND_LIMITED_5);
const InterpolationScheme BL6 = InterpolationScheme(BAND_LIMITED_6);
const InterpolationScheme BL7 = InterpolationScheme(BAND_LIMITED_7);
const InterpolationScheme BL_TO_CELL_0 = InterpolationScheme(BAND_LIMITED_CELL_ZERO);
const InterpolationScheme CBFst = InterpolationScheme(CUBIC_INTERP_FIRST);
const InterpolationScheme CBMid = InterpolationScheme(CUBIC_INTERP_MIDDLE);
const InterpolationScheme CBLst = InterpolationScheme(CUBIC_INTERP_LAST);

/**
 * @brief Determine the appropriate interpolation scheme to use, given the number of datapoints available and the position to which we wish to interpolate to.
 *
 * @param datapts_in_direction The number of datapoints available along this axis
 * @param interpolation_position Interpolation is to be performed to the midpoint of the datapoints indexed with (interpolation_position-1) and (interpolation_position)
 * @param interpolation_methods The preferred interpolation methods to use
 * @return const InterpolationScheme&
 */
const InterpolationScheme &best_scheme(int datapts_in_direction, int interpolation_position,
                                       PreferredInterpolationMethods interpolation_methods =
                                               PreferredInterpolationMethods::BandLimited);
