/**
 * @file openandorder.cpp
 * @brief Launch and file IO
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

#define NMATRICES 49              //< number of input matrices
#define NOUTMATRICES_WRITE 23     //< number of output matrices to be written to output file
#define NOUTMATRICES_WRITE_ALL 25 //< number of output matrices to be written to output file
#define NOUTMATRICES_PASSED 31    //< number of output matrices passed by mexFunction

using namespace std;


int main(int nargs, char *argv[]){

  // Set the logging level with a compile-time define for debugging
  #if SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_DEBUG
    spdlog::set_level(spdlog::level::debug);
  #elif SPDLOG_ACTIVE_LEVEL == SPDLOG_LEVEL_INFO
    spdlog::set_level(spdlog::level::info);
  #endif

  /*
    There are two cases to consider, when the fdtdgrid matrix is specified in a separate mat file
    and when it is in the same file as the other matrices.
  */

  const char *matrixnames[] = {"fdtdgrid","Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
  const char *matrixnames_infile[] = {"Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
  const char *matrixnames_gridfile[] = {"fdtdgrid"};
  const char *outputmatrices_all[] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","vertices","camplitudes","facets","maxresfield","Id","fieldsample","campssample"};
  const char *outputmatrices[] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","camplitudes","maxresfield","Id","fieldsample","campssample"};
  int  matricestosave_all[]   = {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};
  int  matricestosave[]   = {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,23,25,26,27,28};
  mxArray *plhs[NOUTMATRICES_PASSED];
  const mxArray *matrixptrs[NMATRICES];

  auto args = ArgumentParser::parse_args(nargs, argv);
  check_files_can_be_accessed(args);

  //now it is safe to use matlab routines to open the file and order the matrices
  if(!args.has_grid_filename())
    openandorder(args.input_filename(), (char **)matrixnames, matrixptrs, NMATRICES);
  else{
    openandorder(args.input_filename(), (char **)matrixnames_infile, matrixptrs+1, NMATRICES-1);
    openandorder(args.grid_filename(), (char **)matrixnames_gridfile, matrixptrs, 1);
  }

  //now run the time propagation code
  moaf(NOUTMATRICES_PASSED, (mxArray **)plhs, NMATRICES, (const mxArray **)matrixptrs);

  if( !args.have_flag("-m") ){ //prints vertices and facets
    saveoutput(plhs, matricestosave_all, (char **)outputmatrices_all, NOUTMATRICES_WRITE_ALL, args.output_filename());
  }
  else{ // minimise the file size by not printing vertices and facets
    saveoutput(plhs, matricestosave, (char **)outputmatrices, NOUTMATRICES_WRITE, args.output_filename());
  }

  return 0;
}

void openandorder(const char *mat_filename, char **matrix_names, const mxArray **matrix_ptrs, int n_matrices){

  auto expected = MatrixCollection(matrix_names, n_matrices);
  auto actual = MatFileMatrixCollection(mat_filename);

  actual.check_has_at_least_as_many_matrices_as(expected);

  assign_matrix_pointers(expected, actual, matrix_ptrs);

  //fprintf(stderr, "Got all %d matrices\n",nmatrices);
  //fprintf(stderr, "%d\n",mxGetFieldNumber( (mxArray *)matrixptrs[0], "I_tot"));

  if (!are_equal(matrix_names[0], "fdtdgrid")){
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

void saveoutput(mxArray **plhs, const int *matricestosave, char *matrixnames[], int nmatrices, const char *outputfilename){

  auto outfile = matOpen(outputfilename, "w7.3");
  if (outfile == nullptr){
    throw runtime_error("Unable to open file output file" + string(outputfilename));
  }

  //now iterate through the matrices, set names and add to matfile
  for(int i=0; i<nmatrices; i++){
    auto mpv_out = matPutVariable(outfile, matrixnames[i], (mxArray *)plhs[matricestosave[i]]);

    if(mpv_out){
      auto fp = matGetFp(outfile);
      fprintf(stderr, "Could not write array %s to %s (%d,%d,%d)\n",
              matrixnames[i], outputfilename, mpv_out, feof(fp), ferror(fp));
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

      if(are_equal(actual_matrix_name, expected_matrix_name)){
        //fprintf(stderr,"Got %s/%s (%d)\n",mat_matrixnames[j],matrixnames[i],j);
        auto pointer = matGetVariable(actual.mat_file, actual_matrix_name);

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
