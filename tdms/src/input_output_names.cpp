#include "input_output_names.h"

using namespace std;

namespace tdms_matrix_names {
    const vector<string> matrixnames_infile = {"Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
    const vector<string> matrixnames_gridfile = {"fdtdgrid"};
    const vector<string> matrixnames_input_with_grid = {"fdtdgrid","Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
    const vector<string> outputmatrices_all = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","vertices","camplitudes","facets","maxresfield","Id","fieldsample","campssample"};
    const vector<string> outputmatrices = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","camplitudes","maxresfield","Id","fieldsample","campssample"};
    const int matricestosave_all[NOUTMATRICES_WRITE_ALL] = {0,  1,  2,  3,  4,  5,  10, 11, 12,
                                                            13, 14, 15, 16, 17, 18, 19, 20, 21,
                                                            22, 23, 24, 25, 26, 27, 28};
    const int matricestosave[NOUTMATRICES_WRITE] = {0,  1,  2,  3,  4,  5,  10, 11, 12, 13, 14, 15,
                                                    16, 17, 18, 19, 20, 21, 23, 25, 26, 27, 28};
}
