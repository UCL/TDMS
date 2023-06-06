/**
 * @file vector_typedefs.h
 * @author William Graham
 * @brief Typedefs of std::vector for readability improvements across the TDMS
 * codebase.
 */
#pragma once

#include <vector>

/*! Vector of frequencies to extract field & phasors at */
typedef std::vector<double> FrequencyExtractVector;

/* An int array containing the MATLAB indices of field components.

The ints in this array correspond to the underlying values of the
FieldComponents enum: components = [1, 2, 6] corresponds to extracting Ex, Ey,
and Hz at the vertices, for example.
*/
typedef std::vector<int> FieldComponentsVector;
