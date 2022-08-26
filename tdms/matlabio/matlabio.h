#pragma once
#include <algorithm>
#include <complex>
#include <string>
#include "mat_io.h"

/**
 * Casts a 4-dimensional array such that it may be indexed according to the usual array indexing
 * scheme array[l,k,j,i].
 * @param array is a point to a matlab 4 dimensional array
 * @param nrows the number of rows in the array
 * @param ncols the number of columns in the array
 * @param nlayers the number of layers, each of dimension nrows*ncols
 * @param nblocks the number of blocks, each of dimension nrows*ncols*nlayers
 */
double**** cast_matlab_4D_array(double *array, int nrows, int ncols, int nlayers, int nblocks);

void free_cast_matlab_4D_array(double ****castArray, int nlayers, int nblocks);

/**
 * Casts a 3-dimensional array such that it may be indexed according to the usual array indexing
 * scheme array[k,j,i].
 * @param array is a point to a matlab 4 dimensional array
 * @param nrows the number of rows in the array
 * @param ncols the number of columns in the array
 * @param nlayers the number of layers, each of dimension nrows*ncols
 */
template<typename T, typename S>
T*** cast_matlab_3D_array(T *array, S nrows, S ncols, S nlayers){

  T ***p;
  nlayers = std::max(nlayers, 1);
  p = (T ***)malloc((unsigned) (nlayers*sizeof(T **)));

  for(S k =0; k<nlayers; k++){
    p[k] = (T **)malloc((unsigned) (ncols*sizeof(T *)));
  }
  for(S k =0; k<nlayers; k++)
    for(S j =0; j<ncols; j++){
      p[k][j] = (array + k*nrows*ncols+ j*nrows);
    }

  return p;
};

template<typename T>
void free_cast_matlab_3D_array(T ***castArray, int nlayers){
  for(int k =0; k<nlayers;k++)
    free(castArray[k]);
  free(castArray);
}

/**
 * Casts a 2-dimensional array such that it may be indexed according to the usual array indexing
 * scheme array[j,i].
 * @param array is a point to a matlab 4 dimensional array
 * @param nrows the number of rows in the array
 * @param ncols the number of columns in the array
 */
template<typename T, typename S>
T** cast_matlab_2D_array(T *array, S nrows, S ncols){

  T **p;
  p = (T **)malloc((unsigned) (ncols*sizeof(T *)));

  for(S j =0; j<ncols;j++){
    p[j] = (array + j*nrows);
  }
  return p;
};

template<typename T>
void free_cast_matlab_2D_array(T **castArray){
  free(castArray);
}

void assert_is_struct(const mxArray* ptr, const std::string &name);

void assert_num_fields_equals(int num, const mxArray* ptr, const std::string &name);

void assert_is_struct_with_n_fields(const mxArray* ptr, int num, const std::string &name);

/**
 * Get a pointer to a matrix within a struct with a given name. Throws a runtime error if the
 * resulting tensor is not two dimensional.
 * @param ptr Pointer to the struct
 * @param name Name of the attribute
 * @param struct_name Name of the struct
 * @return Pointer to the matrix
 */
mxArray* ptr_to_matrix_in(const mxArray* ptr, const std::string &name, const std::string &struct_name);

/**
 * Get a pointer to an array within a struct with a given name. Throws a runtime error if the
 * resulting tensor is not a vector with dimensions 1xN.
 * @param ptr Pointer to the struct
 * @param name Name of the attribute
 * @param struct_name Name of the struct
 * @return Pointer to the vector
 */
mxArray* ptr_to_vector_in(const mxArray* ptr, const std::string &name, const std::string &struct_name);

mxArray* ptr_to_vector_or_empty_in(const mxArray* ptr, const std::string &name, const std::string &struct_name);

/**
 * Get a double defined in a matlab array given as a pointer
 * @param ptr Pointer to a matlab array
 * @param name Name of the value, for helpful thrown exceptions if the pointer is not to a double
 * @return Value of the double
 */
double double_in(const mxArray* ptr, const std::string &name);

int int_cast_from_double_in(const mxArray* ptr, const std::string &name);

bool bool_cast_from_double_in(const mxArray* ptr, const std::string &name);

/**
 * Get the (C++) string defined in a matlab array given as a pointer
 * @param ptr Pointer to a matlab array
 * @param name Name of the field which this pointer corresponds to
 * @return The string
 */
std::string string_in(const mxArray* ptr, const std::string &name);
