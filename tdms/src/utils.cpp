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

int max(int a, int b, int c) { return std::max(std::max(a, b), c); }
