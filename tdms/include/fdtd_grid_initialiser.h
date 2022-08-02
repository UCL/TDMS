#include "mat_io.h"
#include "vector"


class fdtdGridInitialiser {

private:
    const mxArray *pointer;
    const char *mat_filename;
    std::vector<mwSize> dimensions;

    int value_of_attribute(const std::string& key);

public:
    fdtdGridInitialiser(const mxArray *fdtd_pointer, const char* mat_filename);

    void add_tensor(const std::string &name);
};
