/**
 * A collection of functions for tensor initialisation given matlab pointers
 */
#pragma once
#include <string>
#include <vector>
#include "vector_collection.h"
#include "mat_io.h"
#include "field.h"


/**
 * Initialise the grid tensors/arrays, including the electric and magnetic split fields and the
 * materials array.
 * @param ptr Pointer to the matlab struct
 * @param E_s Electric split field
 * @param H_s Magnetic split field
 * @param materials Materials array
 */
void init_grid_tensors(const mxArray *ptr, SplitField &E_s, SplitField &H_s, uint8_t*** &materials);