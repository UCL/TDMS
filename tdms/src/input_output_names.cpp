#include "input_output_names.h"

namespace tdms_matrix_names {
    const char *matrixnames[NMATRICES] = {"fdtdgrid","Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
    const char *matrixnames_infile[NMATRICES-1] = {"Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
    const char *matrixnames_gridfile[1] = {"fdtdgrid"};
    const char *outputmatrices_all[NOUTMATRICES_WRITE_ALL] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","vertices","camplitudes","facets","maxresfield","Id","fieldsample","campssample"};
    const char *outputmatrices[NOUTMATRICES_WRITE] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","camplitudes","maxresfield","Id","fieldsample","campssample"};
    int matricestosave_all[NOUTMATRICES_WRITE_ALL] = {0,  1,  2,  3,  4,  5,  10, 11, 12,
                                                      13, 14, 15, 16, 17, 18, 19, 20, 21,
                                                      22, 23, 24, 25, 26, 27, 28};
    int matricestosave[NOUTMATRICES_WRITE] = {0,  1,  2,  3,  4,  5,  10, 11, 12, 13, 14, 15,
                                              16, 17, 18, 19, 20, 21, 23, 25, 26, 27, 28};
}
