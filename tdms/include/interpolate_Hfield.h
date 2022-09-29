/**
 * @file interpolate_Hfield.h
 * @brief Interpolation functions for each of the H-field components.
 * 
 * Decided to do with 3 separate functions for readability; as this avoids
 * switches with a lot of cases within functions, and the interpolation scheme
 * might be different in each direction anyway, so there is no point having a
 * centralised function for all of them.
 * 
 * Additionally, generalised functions for performing the 2D interpolation is
 * not pretty code to look at (although is theoretically possible if we wanted a
 * slight speedup).  This forces us to pass in I,J,K to _every_ component
 * function, rather than only the number of cells in the appropriate direction.
 * 
 * If we want to interpolate all the H-field components, there are functions for
 * this which call the individual component functions separately.
 */
#pragma once

void interpolateTimeDomainHx(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx);
void interpolateTimeDomainHy(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx);
void interpolateTimeDomainHz(double ***Hxy, double ***Hxz, int i, int j, int k, int I, int J, int K, double *Hx);
void interpolateTimeDomainHField(double ***Hxy, double ***Hxz, double ***Hyx,
                                 double ***Hyz, double ***Hzx, double ***Hzy,
                                 int i, int j, int k, int I, int J, int K,
                                 double *Hx, double *Hy, double *Hz);

bool determine_second_dim(int n_cells_d0, int n_cells_d1, int cid0, int cid1);