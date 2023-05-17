/**
 * @file input_output_names.h
 * @brief Constant values associated with the format of the input and output
 * files and matrices
 */
#pragma once

#include <string>
#include <vector>

#define NMATRICES 48//< number of input matrices
#define NOUTMATRICES_WRITE                                                     \
  23//< number of output matrices to be written to output file
#define NOUTMATRICES_WRITE_ALL                                                 \
  25//< number of output matrices to be written to output file
#define NOUTMATRICES_PASSED                                                    \
  31//< number of output matrices passed by mexFunction

/*
These are the names of the arrays that we expect to recieve from the input file
(and posibly gridfile) that is passed to tdms on the command line. In the case
where we are given an input file and a gridfile, the array names we expect to
recieve from each is also recorded. In the case where the -m (compressed output)
flag is passed, we do not save all the possible output arrays and so have a
complete list of outputs and a shortened list of those which are saved when
using -m.
*/
namespace tdms_matrix_names {
const std::vector<std::string> matrixnames_infile = {
        "Cmaterial",
        "Dmaterial",
        "C",
        "D",
        "freespace",
        "disp_params",
        "delta",
        "interface",
        "Isource",
        "Jsource",
        "Ksource",
        "grid_labels",
        "omega_an",
        "to_l",
        "hwhm",
        "Dxl",
        "Dxu",
        "Dyl",
        "Dyu",
        "Dzl",
        "Dzu",
        "Nt",
        "dt",
        "tind",
        "sourcemode",
        "runmode",
        "exphasorsvolume",
        "exphasorssurface",
        "intphasorssurface",
        "phasorsurface",
        "phasorinc",
        "dimension",
        "conductive_aux",
        "dispersive_aux",
        "structure",
        "f_ex_vec",
        "exdetintegral",
        "f_vec",
        "Pupil",
        "D_tilde",
        "k_det_obs_global",
        "air_interface",
        "intmatprops",
        "tdfield",
        "tdfdir",
        "fieldsample",
        "campssample"};//< All matrices that we could expect to get from an
                       // input file with a separate gridfile
const std::vector<std::string> matrixnames_input_with_grid = {
        "fdtdgrid",
        "Cmaterial",
        "Dmaterial",
        "C",
        "D",
        "freespace",
        "disp_params",
        "delta",
        "interface",
        "Isource",
        "Jsource",
        "Ksource",
        "grid_labels",
        "omega_an",
        "to_l",
        "hwhm",
        "Dxl",
        "Dxu",
        "Dyl",
        "Dyu",
        "Dzl",
        "Dzu",
        "Nt",
        "dt",
        "tind",
        "sourcemode",
        "runmode",
        "exphasorsvolume",
        "exphasorssurface",
        "intphasorssurface",
        "phasorsurface",
        "phasorinc",
        "dimension",
        "conductive_aux",
        "dispersive_aux",
        "structure",
        "f_ex_vec",
        "exdetintegral",
        "f_vec",
        "Pupil",
        "D_tilde",
        "k_det_obs_global",
        "air_interface",
        "intmatprops",
        "tdfield",
        "tdfdir",
        "fieldsample",
        "campssample"};//< Matrices we expect to get from an input file that
                       // also contains grid information
const std::vector<std::string> matrixnames_gridfile = {
        "fdtdgrid"};//< Matrices we expect to obtain from a separate gridfile
const std::vector<std::string> outputmatrices_all = {
        "Ex_out",     "Ey_out",      "Ez_out",   "Hx_out",
        "Hy_out",     "Hz_out",      "x_out",    "y_out",
        "z_out",      "Ex_i",        "Ey_i",     "Ez_i",
        "Hx_i",       "Hy_i",        "Hz_i",     "x_i",
        "y_i",        "z_i",         "vertices", "camplitudes",
        "facets",     "maxresfield", "Id",       "fieldsample",
        "campssample"};//< All output matrices we might want to write
const std::vector<std::string> outputmatrices = {
        "Ex_out", "Ey_out",      "Ez_out",     "Hx_out",      "Hy_out",
        "Hz_out", "x_out",       "y_out",      "z_out",       "Ex_i",
        "Ey_i",   "Ez_i",        "Hx_i",       "Hy_i",        "Hz_i",
        "x_i",    "y_i",         "z_i",        "camplitudes", "maxresfield",
        "Id",     "fieldsample", "campssample"};//< Output matrices we want to
                                                // write in a compressed (-m)
                                                // output
}// namespace tdms_matrix_names
