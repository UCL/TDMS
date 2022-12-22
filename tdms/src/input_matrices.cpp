#include "input_matrices.h"

#include <string.h>
#include <stdexcept>

#include <spdlog/spdlog.h>

#include "fdtd_grid_initialiser.h"
#include "utils.h"
#include "matlabio.h"

using namespace std;
using namespace tdms_matrix_names;

int InputMatrices::index_from_matrix_name(const string &matrix_name) {
  auto position = find(matrixnames_input_with_grid.begin(), matrixnames_input_with_grid.end(), matrix_name);
  if (position == matrixnames_input_with_grid.end()) {
    // could not find the matrix name in the list of expected input matrices
    throw runtime_error(matrix_name + " not found in matrixnames_input_with_grid");
  }
  return distance(matrixnames_input_with_grid.begin(), position);
}

void InputMatrices::set_from_input_file(const char *mat_filename) {
  MatrixCollection infile_expected(matrixnames_input_with_grid);
  MatFileMatrixCollection infile_contains(mat_filename);
  spdlog::info("Input file: " + string(mat_filename) + " | No gridfile supplied");

  // check that the input file actually has enough matrices for what we're expecting
  infile_contains.check_has_at_least_as_many_matrices_as(infile_expected);
  // assign the pointers to the matrices
  assign_matrix_pointers(infile_expected, infile_contains);

  // Add fields to fdtdgrid
  auto initialiser =
          fdtdGridInitialiser(matrix_pointers[index_from_matrix_name("fdtdgrid")], mat_filename);
  const vector<string> fdtdgrid_element_names = {"Exy", "Exz", "Eyx", "Eyz", "Ezx", "Ezy",
                                                  "Hxy", "Hxz", "Hyx", "Hyz", "Hzx", "Hzy"};

  for (const auto &name : fdtdgrid_element_names) { initialiser.add_tensor(name); }

  // validate the input arguments
  validate_assigned_pointers();
}
void InputMatrices::set_from_input_file(const char *mat_filename, const char *gridfile) {
  MatrixCollection infile_expected(matrixnames_infile);
  MatFileMatrixCollection infile_contains(mat_filename);
  spdlog::info("Input file: " + string(mat_filename) + " | Gridfile: " + string(gridfile));

  // first extract fdtdgrid from the gridfile provided
  MatrixCollection gridfile_expected(matrixnames_gridfile);
  MatFileMatrixCollection gridfile_contains(gridfile);
  gridfile_contains.check_has_at_least_as_many_matrices_as(gridfile_expected);

  // this is the matrix name we are searching for in the gridfile
  string fdtd_matrix_search_string = gridfile_expected.matrix_names[0];
  // check that the gridfile actually contains the fdtdgrid
  auto position = find(gridfile_contains.matrix_names.begin(), gridfile_contains.matrix_names.end(),
                       fdtd_matrix_search_string);
  if (position == matrixnames_input_with_grid.end()) {
    throw runtime_error(fdtd_matrix_search_string + " not found in gridfile");
  }
  // attempt to set the pointer
  auto pointer = matGetVariable(gridfile_contains.mat_file, fdtd_matrix_search_string.c_str());
  if (pointer == nullptr) {
    throw runtime_error("Could not get pointer to " + fdtd_matrix_search_string);
  }
  set_matrix_pointer("fdtdgrid", pointer);

  // now pick up the other matrix inputs from the usual input file
  // check that the input file actually has enough matrices for what we're expecting
  infile_contains.check_has_at_least_as_many_matrices_as(infile_expected);
  // assign the pointers to the matrices
  assign_matrix_pointers(infile_expected, infile_contains);

  // Add fields to fdtdgrid
  auto initialiser =
          fdtdGridInitialiser(matrix_pointers[index_from_matrix_name("fdtdgrid")], mat_filename);
  const vector<string> fdtdgrid_element_names = {"Exy", "Exz", "Eyx", "Eyz", "Ezx", "Ezy",
                                                  "Hxy", "Hxz", "Hyx", "Hyz", "Hzx", "Hzy"};

  for (const auto &name : fdtdgrid_element_names) { initialiser.add_tensor(name); }

  // validate the input arguments
  validate_assigned_pointers();
}

void InputMatrices::assign_matrix_pointers(MatrixCollection &expected,
                                           MatFileMatrixCollection &actual) {
  // for each matrix we expect to recieve
  for(string &expected_matrix : expected.matrix_names) {
    // look for this matrix name in the input file, throw an error if it is not present
    auto position = find(actual.matrix_names.begin(), actual.matrix_names.end(), expected_matrix);
    if (position == actual.matrix_names.end()) {
      throw runtime_error(expected_matrix + " not found in input file");
    }
    // if it's present, attempt to set the corresponding pointer
    auto pointer = matGetVariable(actual.mat_file, expected_matrix.c_str());
    if (pointer == nullptr) {
      throw runtime_error("Could not get pointer to " + expected_matrix);
    } else {
      set_matrix_pointer(expected_matrix, pointer);
    }
  }
}

void InputMatrices::validate_assigned_pointers() {
  // certain arrays must be structure arrays
  assert_is_struct(matrix_pointers[index_from_matrix_name("fdtdgrid")], "fdtdgrid");
  assert_is_struct(matrix_pointers[index_from_matrix_name("Cmaterial")], "Cmaterial");
  assert_is_struct(matrix_pointers[index_from_matrix_name("Dmaterial")], "Dmaterial");
  assert_is_struct(matrix_pointers[index_from_matrix_name("C")], "C");
  assert_is_struct(matrix_pointers[index_from_matrix_name("D")], "D");

  // other arrays must contain a specific number of fieldnames
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("freespace")], 6,
                                  "freespace");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("disp_params")], 3,
                                  "disp_params");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("delta")], 3, "delta");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("interface")], 6,
                                  "interface");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("grid_labels")], 3,
                                  "grid_labels");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("conductive_aux")], 3,
                                  "conductive_aux");
}
