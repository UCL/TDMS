/**
 * @file openandorder.cpp
 * @brief Passing arguments and file IO
 *
 * Code for processing command line arguments, opening input files,  passing
 * matrices to the mexFunction and writing the output to the specified output
 * file.
 */
#include "openandorder.h"

#include <cstdio>
#include <stdexcept>

#include <spdlog/spdlog.h>

#include "utils.h"
#include "fdtd_grid_initialiser.h"
#include "iterator.h"

using namespace std;

void saveoutput(mxArray **plhs, const int *_matricestosave,
                const std::vector<std::string> &matrixnames, int nmatrices,
                const char *outputfilename) {

  auto outfile = matOpen(outputfilename, "w7.3");
  if (outfile == nullptr){
    throw runtime_error("Unable to open file output file" + string(outputfilename));
  }

  //now iterate through the matrices, set names and add to matfile
  for(int i=0; i<nmatrices; i++){
    auto mpv_out = matPutVariable(outfile, matrixnames[i].c_str(), (mxArray *)plhs[_matricestosave[i]]);

    if(mpv_out){
      auto fp = matGetFp(outfile);
      spdlog::info("Could not write array {0:s} to {1:s} ({2:d},{3:d},{4:d})", matrixnames[i].c_str(), outputfilename, mpv_out, feof(fp), ferror(fp));
    }
  }

  matClose(outfile);
}

void check_files_can_be_accessed(ArgumentNamespace &args){

  for (const auto& filename: args.input_filenames()){
    assert_can_open_file(filename.c_str(), "r");
  }

  assert_can_open_file(args.output_filename(), "a+");
}
