/*****************************************************************
 *  Application.:  matlab data structure manipulation
 *  Description.:  code for processing the matlab data structures
 ******************************************************************/
#include "matlabio.h"


using namespace std;

double****cast_matlab_4D_array(double *array, int nrows, int ncols, int nlayers, int nblocks){
  double ****p;
  int j,k,l;

  if( nlayers==0 )
	  nlayers++;
  
  if( nblocks==0 )
	  nblocks++;
  
  p = (double ****)malloc((unsigned) (nblocks*sizeof(double ***)));
  for(l=0; l<nblocks; l++)
    p[l] = (double ***)malloc((unsigned) (nlayers*sizeof(double **)));
  
  for(l=0; l<nblocks; l++)
    for(k =0; k<nlayers;k++)
      p[l][k] = (double **)malloc((unsigned) (ncols*sizeof(double *)));

  for(l=0; l<nblocks; l++)
    for(k =0; k<nlayers;k++)
      for(j =0; j<ncols;j++)
	p[l][k][j] = (array + l*nrows*ncols*nlayers + k*nrows*ncols + j*nrows);
  return p;

}

void free_cast_matlab_4D_array(double ****castArray, int nlayers, int nblocks){
  for(int l =0; l<nblocks; l++)
    for(int k =0; k<nlayers;k++)
      free(castArray[l][k]);

  for(int l =0; l<nblocks; l++)
    free(castArray[l]);
  free(castArray);
}

void assert_is_struct(const mxArray* ptr, const std::string &name){
  if (!mxIsStruct(ptr)) {
    throw runtime_error("Pointer " + name + " was expected to be a structure but was not");
  }
}

void assert_num_fields_equals(int num, const mxArray* ptr, const std::string &name){

  auto num_fields = mxGetNumberOfFields(ptr);

  if (num_fields != num) {
    throw runtime_error(name + " should have " + to_string(num) + " members, it has "
                        + to_string(num_fields));
  }
}

void assert_is_struct_with_n_fields(const mxArray* ptr, int num, const std::string &name){
  assert_is_struct(ptr, name);
  assert_num_fields_equals(num, ptr, name);
}

mxArray* ptr_to_nd_array_in(const mxArray* ptr, int n, const std::string &name, const std::string &struct_name){
  auto element = mxGetField((mxArray *) ptr, 0, name.c_str());

  if (mxGetNumberOfDimensions(element) != n) {
    throw runtime_error("Incorrect dimension on " + struct_name + "." + name + ". Required " +
                        to_string(n) + "D");
  }
  return element;
}

mxArray* ptr_to_matrix_in(const mxArray* ptr, const string &name, const string &struct_name){
  return ptr_to_nd_array_in(ptr, 2, name, struct_name);
}

mxArray* ptr_to_vector_in(const mxArray* ptr, const string &name, const string &struct_name){

  auto element = ptr_to_matrix_in(ptr, name, struct_name);

  if (mxGetDimensions(element)[0] != 1) {
    throw runtime_error("Incorrect dimension " + struct_name + "." + name + ". Required 1D");
  }
  return element;
}

mxArray* ptr_to_vector_or_empty_in(const mxArray* ptr, const string &name, const string &struct_name){

  auto element = ptr_to_matrix_in(ptr, name, struct_name);
  auto n = mxGetDimensions(element)[0];

  if (!(n == 1 || n == 0)) {
    throw runtime_error("Incorrect dimension on " + struct_name + "." + name + ". Required 1D or 0D");
  }
  return element;
}

double double_in(const mxArray* ptr, const string &name){

  if (mxIsDouble(ptr) && mxGetNumberOfElements(ptr) == 1) {
    return *mxGetPr(ptr);
  }
  throw runtime_error(name + " was expected to be a double but was not");
}

int int_cast_from_double_in(const mxArray* ptr, const std::string &name){
  return (int) double_in(ptr, name);
}

bool bool_cast_from_double_in(const mxArray* ptr, const std::string &name){
  return (bool) double_in(ptr, name);
}

string string_in(const mxArray *ptr, const string &name) {

  if (mxIsChar(ptr)){
    auto n = 1 + (int) mxGetNumberOfElements(ptr);
    auto c_str = (char *) malloc(n * sizeof(char));
    mxGetString(ptr, c_str, n);
    return c_str;
  }
  throw runtime_error(name + " was expected to be a string but was not");
}
