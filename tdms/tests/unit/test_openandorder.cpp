#include <catch2/catch_test_macros.hpp>
#include "argument_parser.h"

using namespace std;


TEST_CASE("Test argument parsing") {

  const char *input_args[] = {"tdms", "-h"};
  auto args = ArgumentNamespace(2, const_cast<char **>(input_args));
  REQUIRE(args.have_flag("-h"));
}
