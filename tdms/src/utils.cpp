#include "utils.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

using namespace std;

void assert_can_open_file(const char *filename, const char *mode) {

  auto file = fopen(filename, mode);

  if (file == nullptr) {
    throw runtime_error("Unable to open file " + string(filename));
  }

  fclose(file);
}

bool are_equal(const char *a, const char *b) { return strcmp(a, b) == 0; }

vector<int>
tdms_vector_utils::to_vector_int(const vector<double> &ints_that_are_doubles) {
  vector<int> values_as_ints;
  values_as_ints.resize(ints_that_are_doubles.size());
  for (unsigned int i = 0; i < ints_that_are_doubles.size(); i++) {
    values_as_ints[i] = static_cast<int>(ints_that_are_doubles[i]);
  }
  return values_as_ints;
}
