/**
 * @file numerical_derivative.h
 * @brief Functions to calculate the numerical derivatives.
 *
 * This code is the main dependency on FFTW3.
 */
#pragma once
#include <fftw3.h>

/**
 * @brief Multiply two arrays of complex numbers element-wise.
 *
 * Performs element-wise complex multiplication of the array a with the array b.
 * a and b must be must be of equal length.
 * Stores the result in c.
 * Internally used by first_derivative.
 *
 * @param[in] a Array of complex numbers to multiply with those in b.
 * @param[in] b Array of complex numbers to multiply with those in a.
 * @param[out] c The results of the element-wise multiplications a .* b.
 * @param[in] len The length of the arrays.
 */
void complex_mult_vec(fftw_complex *a, fftw_complex *b, fftw_complex *c,
                      int len);

/**
 * @brief Initialise the coefficients required to simultaneously perform
 * differentiation and shifting by amount delta, using a forward and backward
 * FFT.
 *
 * @param[in] delta The fraction of the spatial step.
 * @param[out] Dk Buffer to write the coefficients to.
 * @param[in] N The number of elements in Dk.
 */
void init_diff_shift_op(double delta, fftw_complex *Dk, int N);

/**
 * @brief Calculate the first derivative of a sampled function.
 * @note in_pb_pf must be the buffer which is the input for both plans pf and pb.
 * Likewise, out_pb_pf must be the output for both plans pf and pb.
 *
 * @param[in] in_pb_pf The buffer containing the data to be differentiated.
 * 	    @warning This buffer will be overwritten as part of the computation.
 * @param[out] out_pb_pf The buffer which will contain the computed derivative.
 * @param[inout] Dk Buffer to write the coefficients for performing differentiation and shifting in Fourier space.
 * @param[in] N Number of elements in buffers.
 * @param[in] pf The plan for forward FFT.
 * @param[in] pb The plan for backward FFT.
 */
void first_derivative(fftw_complex *in_pb_pf, fftw_complex *out_pb_pf,
                      fftw_complex *Dk, int N, fftw_plan pf, fftw_plan pb);
