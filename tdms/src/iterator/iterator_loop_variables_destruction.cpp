/**
 * @file iterator_loop_variables_destructor.cpp
 * @brief Destructor method for Iterator_LoopVariables class, and the methods called in the destructor
 */
#include "iterator_loop_variables.h"

#include <spdlog/spdlog.h>

void Iterator_LoopVariables::free_Id_memory() {
  free_cast_matlab_2D_array(Idx_re);
  free_cast_matlab_2D_array(Idx_im);
  free_cast_matlab_2D_array(Idy_re);
  free_cast_matlab_2D_array(Idy_im);
  for (int ifx = 0; ifx < f_ex_vec.size(); ifx++) {
    free(Idx[ifx]);
    free(Idy[ifx]);
  }
  free(Idx);
  free(Idy);
}

void Iterator_LoopVariables::free_iwave_memory() {
  free_cast_matlab_2D_array(iwave_lEx_Rbs);
  free_cast_matlab_2D_array(iwave_lEx_Ibs);
  free_cast_matlab_2D_array(iwave_lEy_Rbs);
  free_cast_matlab_2D_array(iwave_lEy_Ibs);
  free_cast_matlab_2D_array(iwave_lHx_Rbs);
  free_cast_matlab_2D_array(iwave_lHx_Ibs);
  free_cast_matlab_2D_array(iwave_lHy_Rbs);
  free_cast_matlab_2D_array(iwave_lHy_Ibs);
}

Iterator_LoopVariables::~Iterator_LoopVariables() {
  /* ORIGINAL TEAR-DOWN ORDER WAS:
  - mxFree(mx_surface_vertices) (if params.exphasorssurface && complete run mode)
  - ~SurfacePhasors (if params.exphasorssurface && complete run mode)
  - ~VertexPhasors
  - Id-related variable cleanup (free_Id_memory)
  - {I,J,K}-source tear-down (happens in ~Iterator_ObjectsFromInfile now)
  - iwave variables (always cleared)
  - free_cast materials (~Iterator_IndependentObjectsFromInfile handles this)
  - free PSTD memory allocation
  - free {E,H}-norm, and dims and label_dims. The later two are no longer assigned by malloc so we don't need it here
  - free E_copy_data_placeholders (if steadystate and complete run mode)

  This is now slightly out of order with how C++ will call superclass destructors after this one, but by construction none of the superclasses should ever be touching the same memory, whilst simultaneously taking care of their own.
  */
  spdlog::info("Destroying Iterator_LoopVariables");

  if (params.exphasorssurface && params.run_mode == RunMode::complete) {
    mxFree(mx_surface_vertices);
    // ~SurfacePhasors object will clean up the remaining surface-phasor memory
  }
  // ~VertexPhasors will clean up vertex phasor memory

  // free C++ memory that was cast to the Id output, but is no longer needed
  if (params.exdetintegral) {
    // We have cast the Id variables to a MATLAB array, and must free the memory
    free_Id_memory();
  }

  // {I,J,K}-source are now torn down by ~Iterator_ObjectsFromInfile

  // free the iwave variables (under all circumstances)
  free_iwave_memory();

  // materials torn down by ~Iterator_IndependentObjectsFromInfile

  // free PSTD-unique memory, if we used the PSTD method
  if (solver_method == SolverMethod::PseudoSpectral) {
    // if we used the psuedo-spectral method, we need to delete the malloc'd fftw space for the derivative-shift operators
    fftw_free(dk_e_x);
    fftw_free(dk_e_y);
    fftw_free(dk_e_z);
    fftw_free(dk_h_x);
    fftw_free(dk_h_y);
    fftw_free(dk_h_z);
  }

  // free {E,H}_norm malloc'd memory
  free(E_norm);
  free(H_norm);

  // free E_copy_data_placeholders, which we were using to record the previous phasor state to check convergence,  if we were running a steady-state simulation
  if (params.source_mode == SourceMode::steadystate && params.run_mode == RunMode::complete) {
    mxDestroyArray(E_copy_data_placeholders[0]);
    mxDestroyArray(E_copy_data_placeholders[1]);
    mxDestroyArray(E_copy_data_placeholders[2]);
  }
}
