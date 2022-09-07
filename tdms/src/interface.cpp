#include <string>
#include <algorithm>
#include "matlabio.h"
#include "interface.h"

using namespace std;

InterfaceComponent::InterfaceComponent(const mxArray *ptr, const string &name){

  auto element = ptr_to_matrix_in(ptr, name, "interface");
  auto dims = mxGetDimensions(element);

  if (!(dims[0] == 1 && dims[1] == 2)) {
    throw runtime_error("Incorrect dimension on interface." + name + " (" +
                        to_string((int) dims[0]) + "," +
                        to_string((int) dims[1]) + "," +
                        to_string((int) mxGetNumberOfElements(element)) + ")\n");
  }

  auto array = mxGetPr((mxArray *) element);
  index = max((int)array[0] - 1, 0);   // convert matlab to C indexing, ensuring a min of at least 0
  apply = (bool)array[1];
}
