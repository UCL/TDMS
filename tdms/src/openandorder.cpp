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

#include "utils.h"
#include "fdtd_grid_initialiser.h"
#include "iterator.h"

using namespace std;

void openandorder(const char *mat_filename,
                  const std::vector<std::string>& matrix_names,
                  const mxArray **matrix_ptrs, int n_matrices){

  auto expected = MatrixCollection(matrix_names);
  auto actual = MatFileMatrixCollection(mat_filename);

  actual.check_has_at_least_as_many_matrices_as(expected);

  assign_matrix_pointers(expected, actual, matrix_ptrs);

  //fprintf(stderr, "Got all %d matrices\n",nmatrices);
  //fprintf(stderr, "%d\n",mxGetFieldNumber( (mxArray *)matrixptrs[0], "I_tot"));

  if (matrix_names.back() != "fdtdgrid"){
    return;
  }

  // Add fields to fdtdgrid
  auto initialiser = fdtdGridInitialiser(matrix_ptrs[0], mat_filename);
  const vector<string> fdtdgrid_element_names = {
          "Exy", "Exz", "Eyx", "Eyz", "Ezx", "Ezy",
          "Hxy", "Hxz", "Hyx", "Hyz", "Hzx", "Hzy"
  };

  for(const auto& name : fdtdgrid_element_names){
    initialiser.add_tensor(name);
  }
}

void saveoutput(mxArray **plhs, const int *matricestosave,
                const std::vector<std::string>& matrixnames, int nmatrices,
                const char *outputfilename){

  auto outfile = matOpen(outputfilename, "w7.3");
  if (outfile == nullptr){
    throw runtime_error("Unable to open file output file" + string(outputfilename));
  }

  //now iterate through the matrices, set names and add to matfile
  for(int i=0; i<nmatrices; i++){
    auto mpv_out = matPutVariable(outfile, matrixnames[i].c_str(), (mxArray *)plhs[matricestosave[i]]);

    if(mpv_out){
      auto fp = matGetFp(outfile);
      fprintf(stderr, "Could not write array %s to %s (%d,%d,%d)\n",
              matrixnames[i].c_str(), outputfilename, mpv_out, feof(fp), ferror(fp));
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

void assign_matrix_pointers(MatrixCollection &expected, MatFileMatrixCollection &actual, const mxArray **pointers){

  for(int i=0; i < expected.n_matrices; i++){
    auto expected_matrix_name = expected.matrix_names[i];

    for(int j=0; j < actual.n_matrices; j++){
      auto actual_matrix_name = actual.matrix_names[j];

      if(actual_matrix_name == expected_matrix_name){
        //fprintf(stderr,"Got %s/%s (%d)\n",mat_matrixnames[j],matrixnames[i],j);
        auto pointer = matGetVariable(actual.mat_file, actual_matrix_name.c_str());

        if(pointer == nullptr){
          throw runtime_error("Could not get pointer to "+string(actual_matrix_name));
        }

        pointers[i] = pointer;
        break;
      }
      else if(j==(actual.n_matrices - 1)){  //matrix pointer NOT found
        throw runtime_error("Couldn't find matrix "+string(expected_matrix_name));
      }
    }// j
  }// i
}
