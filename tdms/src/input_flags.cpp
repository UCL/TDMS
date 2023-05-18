#include "input_flags.h"

bool InputFlags::fetch_flag_value(const std::string flag_name,
                                  bool fail_on_not_found) const {
  MATFile *input_file = matOpen(input_filename.c_str(), "r");
  bool return_value;

  mxArray *ptr_to_flag = matGetVariable(input_file, flag_name.c_str());
  if (ptr_to_flag == nullptr) {
    // This flag was not present in the input file.
    // We either return the default value (0, flag not present), or we throw
    // an error if fail_on_not_found is set to true.
    if (fail_on_not_found) {
      throw std::runtime_error(flag_name + " was not present in " +
                               input_filename);
    } else {
      return_value = false;
    }
  } else if (mxIsLogicalScalar(ptr_to_flag)) {
    // The flag is present in the input file, and is a scalar boolean.
    // Return its value
    return_value = (bool) *mxGetPr(ptr_to_flag);
  } else {
    // The flag appears to be present, but is not scalar-valued. Throw
    // exception.
    throw std::runtime_error(flag_name + " is present in " + input_filename +
                             ", but is not scalar.");
  }

  // Cleanup MATLAB structs that are still in memory
  matClose(input_file);
  mxDestroyArray(ptr_to_flag);
  return return_value;
}
