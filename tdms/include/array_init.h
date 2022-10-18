/**
 * @file array_init.h
 * @brief A collection of functions for array initialisation given MATLAB
 * pointers.
 */
#pragma once

#include "arrays.h"
#include "field.h"
#include "mat_io.h"
#include <string>
#include <vector>

/**
 * Initialise the grid tensors/arrays, including the electric and magnetic split fields and the
 * materials array.
 * @param ptr Pointer to the matlab struct
 * @param E_s Electric split field
 * @param H_s Magnetic split field
 * @param materials Materials array
 */
void init_grid_arrays(const mxArray *ptr, SplitField &E_s, SplitField &H_s, uint8_t*** &materials);
