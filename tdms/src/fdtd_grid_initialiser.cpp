#include <string>
#include <stdexcept>
#include "fdtd_grid_initialiser.h"

using namespace std;


fdtdGridInitialiser::fdtdGridInitialiser(const mxArray *fdtd_pointer, const char* fdtd_filename) {
  pointer = fdtd_pointer;
  mat_filename = fdtd_filename;

  dimensions = {value_of_attribute("I_tot") + 1,
                value_of_attribute("J_tot") + 1,
                value_of_attribute("K_tot") + 1};
}


/**
 * Get a value from a integer attribute of fdtdgrid defined in a .mat file
 */
mwSize fdtdGridInitialiser::value_of_attribute(const string& key){

  if( mxGetFieldNumber( (mxArray *)pointer, key.c_str()) == -1 ){
    throw runtime_error(string(mat_filename)+" missing field fdtdgrid."+key);
  }
  auto element = mxGetField( (mxArray *)pointer, 0, key.c_str());

  if (element == nullptr){
    throw runtime_error("Failed to find "+key+" in fdtdgrid");
  }

  auto value = (mwSize) (*mxGetPr(element));
  mxRemoveField((mxArray *)pointer, mxGetFieldNumber(pointer, key.c_str()));

  return value;
}

/**
 * Set a fdtdgrid attribute to a tensor full of zeros
 */
void fdtdGridInitialiser::add_tensor(const string &name){

  mxAddField((mxArray *)pointer, name.c_str());

  auto element = mxGetField((mxArray *)pointer, 0, name.c_str());
  element = mxCreateNumericArray((const mwSize)dimensions.size(),
                                 (const mwSize *)&dimensions.front(),
                                 mxDOUBLE_CLASS, mxREAL);

  mxSetField(( mxArray *)pointer, 0, name.c_str(), element);
}

