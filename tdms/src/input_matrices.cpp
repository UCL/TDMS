#include "input_matrices.h"

#include <string.h>
#include <stdexcept>

#include <spdlog/spdlog.h>

#include "fdtd_grid_initialiser.h"
#include "utils.h"
#include "matlabio.h"

using namespace std;
using namespace tdms_matrix_names;

int index_from_matrix_name(const char *matrix_name) {
    for (int i = 0; i < NMATRICES; i++) {
        if(are_equal(matrix_name, matrixnames[i])) { return i; }
    }
    throw runtime_error("Matrix name " + string(matrix_name) + " not found in matrixnames.");
}

void InputMatrices::set_from_input_file(const char *mat_filename) {
  MatrixCollection infile_expected((char **) matrixnames, NMATRICES);
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
  MatrixCollection infile_expected((char **) matrixnames_infile, NMATRICES - 1);
  MatFileMatrixCollection infile_contains(mat_filename);
  spdlog::info("Input file: " + string(mat_filename) + " | Gridfile: " + string(gridfile));

  // first extract fdtdgrid from the gridfile provided
  MatrixCollection gridfile_expected((char **) matrixnames_gridfile, 1);
  MatFileMatrixCollection gridfile_contains(gridfile);
  gridfile_contains.check_has_at_least_as_many_matrices_as(gridfile_expected);
  // this is the matrix name we are searching for in the gridfile
  char *fdtd_matrix_search_string = gridfile_expected.matrix_names[0];

  // assign the fdtdgrid array
  for (int j = 0; j < gridfile_contains.n_matrices; j++) {
    // the gridfile might have multiple arrays saved in it, but we only want the fdtdgrid
    auto current_matrix_name = gridfile_contains.matrix_names[j];

    if (are_equal(current_matrix_name, fdtd_matrix_search_string)) {
      auto pointer = matGetVariable(gridfile_contains.mat_file, current_matrix_name);

      if (pointer == nullptr) {
        throw runtime_error("Could not get pointer to " + string(current_matrix_name));
      }
      set_matrix_pointer("fdtdgrid", pointer);
      break;
    } else if (j == (gridfile_contains.n_matrices - 1)) {
      // fdtdgrid matrix pointer NOT found
      throw runtime_error("Couldn't find matrix " + string(current_matrix_name));
    }
  }

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
    for (int i = 0; i < expected.n_matrices; i++) {
        auto expected_matrix_name = expected.matrix_names[i];
        // look for this name in the input file
        for (int j = 0; j < actual.n_matrices; j++) {
          auto actual_matrix_name = actual.matrix_names[j];

          if (are_equal(actual_matrix_name, expected_matrix_name)) {
            // if you find the name, get the pointer to the corresponding array
            auto pointer = matGetVariable(actual.mat_file, actual_matrix_name);

            if (pointer == nullptr) {
              // error if we recieve a bad pointer
              throw runtime_error("Could not get pointer to " + string(actual_matrix_name));
            }
            // assign the pointer to the list of pointers
            set_matrix_pointer(actual_matrix_name, pointer);
            break;
          } else if (j == (actual.n_matrices - 1)) {
            //matrix pointer NOT found
            throw runtime_error("Couldn't find matrix " + string(expected_matrix_name));
          }
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
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("disp_params")], 6,
                                  "disp_params");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("delta")], 3, "delta");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("interface")], 6,
                                  "interface");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("grid_labels")], 3,
                                  "grid_labels");
  assert_is_struct_with_n_fields(matrix_pointers[index_from_matrix_name("conductive_aux")], 3,
                                  "conductive_aux");
}
