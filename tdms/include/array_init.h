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


void init_grid_arrays(const mxArray *ptr, SplitField &E_s, SplitField &H_s, uint8_t*** &materials);