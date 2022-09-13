#include <string>
#include <stdexcept>
#include "mat_io.h"
#include "matrix_collection.h"


using namespace std;


MatFileMatrixCollection::MatFileMatrixCollection(const char *filename){

  mat_file = matOpen(filename, "r");
  if(mat_file == nullptr){
    throw runtime_error("Error opening "+string(filename));
  }

  int tmp_n;
  matrix_names = matGetDir(mat_file, &tmp_n);
  n_matrices = tmp_n;
}

MatrixCollection::MatrixCollection(char **names, int number) {

  matrix_names = names;
  n_matrices = number;
}

void MatrixCollection::check_has_at_least_as_many_matrices_as(MatrixCollection &other) {

  if(n_matrices < other.n_matrices){
    throw runtime_error("Not enough matrices in mat file: " +
                        to_string(n_matrices) + " " +
                        to_string(other.n_matrices));
  }
}

MatFileMatrixCollection::~MatFileMatrixCollection() {

  if (mat_file != nullptr){
    matClose(mat_file);
  }
  mxFree(matrix_names);
}
