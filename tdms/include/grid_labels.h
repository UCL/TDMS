#pragma once
#include "matlabio.h"

/**
 * Grid labels hold the cartesian labels of Yee cell, in the x, y and z directions
 */
class GridLabels{
public:
    double *x = nullptr;  // Start of the labels in the x direction
    double *y = nullptr;  //                            y
    double *z = nullptr;  //                            z

    GridLabels() = default;

    explicit GridLabels(const mxArray *ptr);
};
