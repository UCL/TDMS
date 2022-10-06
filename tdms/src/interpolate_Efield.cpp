#include "interpolate_Efield.h"

#include "interpolation_methods.h"

using namespace std;
/*
template<typename T>
void interpolateSplitFieldEx(T ***Exy, T ***Exz, int i, int j, int k, int nI, T *Ex) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nI, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Exy[k][j][i - scheme.number_of_datapoints_to_left + ind] + Exz[k][j][i - scheme.number_of_datapoints_to_left + ind];
    }

    // now run the interpolation scheme and place the result into the output
    *Ex = scheme.interpolate(interp_data);
}
template<typename T>
void interpolateEx(T ***Ex, int i, int j, int k, int nI, T *Ex_interp) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nI, i);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    T interp_data[8];

    // now fill the interpolation data
    // i - (scheme.number_of_datapoints_to_left) is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ex[k][j][i - scheme.number_of_datapoints_to_left + ind];
    }
    // now run the interpolation scheme and place the result into the output
    *Ex_interp = scheme.interpolate(interp_data);
}

void interpolateSplitFieldEy(double ***Eyx, double ***Eyz, int i, int j, int k, int nJ, double *Ey) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nJ, j);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Eyx[k][j - scheme.number_of_datapoints_to_left + ind][i] + Eyz[k][j - scheme.number_of_datapoints_to_left + ind][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ey = scheme.interpolate(interp_data);
}
void interpolateEy(complex<double> ***Ey, int i, int j, int k, int nJ, complex<double> *Ey_interp) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nJ, j);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    complex<double> interp_data[8];

    // now fill the interpolation data
    // j - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ey[k][j - scheme.number_of_datapoints_to_left + ind][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ey_interp = scheme.interpolate(interp_data);
}

void interpolateSplitFieldEz(double ***Ezx, double ***Ezy, int i, int j, int k, int nK, double *Ez) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nK, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    double interp_data[8];

    // now fill the interpolation data
    // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ezx[k - scheme.number_of_datapoints_to_left + ind][j][i] + Ezy[k - scheme.number_of_datapoints_to_left + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ez = scheme.interpolate(interp_data);
}
void interpolateEz(complex<double> ***Ez, int i, int j, int k, int nK, complex<double> *Ez_interp) {
    // determine the interpolation scheme to use
    const interpScheme &scheme = best_interp_scheme(nK, k);
    // prepare input data - if using a cubic scheme we have reserved more memory than necessary but nevermind
    complex<double> interp_data[8];

    // now fill the interpolation data
    // k - scheme.number_of_datapoints_to_left is the index of the Yee cell that plays the role of v0 in the interpolation
    for (int ind = scheme.first_nonzero_coeff; ind <= scheme.last_nonzero_coeff; ind++) {
        interp_data[ind] = Ez[k - scheme.number_of_datapoints_to_left + ind][j][i];
    }

    // now run the interpolation scheme and place the result into the output
    *Ez_interp = scheme.interpolate(interp_data);
}

void interpolateSplitFieldE(double ***Exy, double ***Exz, double ***Eyx,
                            double ***Eyz, double ***Ezx, double ***Ezy,
                            int i, int j, int k, int nI, int nJ, int nK,
                            double *Ex, double *Ey, double *Ez) {
    interpolateSplitFieldEx(Exy, Exz, i, j, k, nI, Ex);
    interpolateSplitFieldEy(Eyx, Eyz, i, j, k, nJ, Ey);
    interpolateSplitFieldEz(Ezx, Ezy, i, j, k, nK, Ez);
}
void interpolateE(complex<double> ***Ex, complex<double> ***Ey, complex<double> ***Ez,
                  int i, int j, int k, int nI, int nJ, int nK,
                  complex<double> *Ex_interp, complex<double> *Ey_interp, complex<double> *Ez_interp) {
    interpolateEx(Ex, i, j, k, nI, Ex_interp);
    interpolateEy(Ey, i, j, k, nJ, Ey_interp);
    interpolateEz(Ez, i, j, k, nK, Ez_interp);
}
*/