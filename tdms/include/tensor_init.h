/**
 * A collection of functions for tensor initialisation given matlab pointers
 */
#include "mat_io.h"
#include "field.h"


void init_grid_tensors(const mxArray *ptr, SplitField &E_s, SplitField &H_s, uint8_t*** &materials);
