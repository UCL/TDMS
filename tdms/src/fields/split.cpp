#include <cstdlib>
#include "field.h"


#include <iostream>

using namespace std;


void SplitFieldComponent::initialise_fftw_plan(int n_threads, int size, EHVec &eh_vec) {

  this->n_threads = n_threads;
  plan_f = (fftw_plan *) malloc(sizeof(fftw_plan *) * n_threads);
  plan_b = (fftw_plan *) malloc(sizeof(fftw_plan *) * n_threads);

  for (int i = 0; i < n_threads; i++){
    plan_f[i] = fftw_plan_dft_1d(size, eh_vec[i], eh_vec[i], FFTW_FORWARD, FFTW_MEASURE);
    plan_b[i] = fftw_plan_dft_1d(size, eh_vec[i], eh_vec[i], FFTW_BACKWARD, FFTW_MEASURE);
  }
}

void SplitFieldComponent::initialise_from_matlab(double ***tensor, Dimensions &dims) {
  this->tensor = tensor;
  this->n_layers = dims[2];
  this->n_cols = dims[1];
  this->n_rows = dims[0];
  this->is_matlab_initialised = true;
}

SplitFieldComponent::~SplitFieldComponent() {

  for (auto plan : {plan_f, plan_b}){
    if (plan == nullptr) continue;

    for (int i = 0; i < n_threads; i++){
      fftw_destroy_plan(plan[i]);
    }
    free(plan);
  }
  // superclass Tensor3D destructor removes tensor memory
}

void SplitField::allocate() {

  for (auto component : {&xy, &xz, &yx, &yz, &zx, &zy}){
    component->allocate(K_tot + 1, J_tot + 1, I_tot + 1);
  }
}

SplitField::SplitField(int I_total, int J_total, int K_total) {

  I_tot = I_total;
  J_tot = J_total;
  K_tot = K_total;
}

void SplitField::zero() {

  for (auto component : {&xy, &xz, &yx, &yz, &zx, &zy}){
    component->zero();
  }
}

void SplitField::initialise_fftw_plan(int n_threads, EHVec &eh_vec) {

  // N_e_y
  int n_x = I_tot + delta_n() + 1;
  int n_y = J_tot + delta_n() + 1;
  int n_z = K_tot + delta_n() + 1;

  xy.initialise_fftw_plan(n_threads, n_y, eh_vec);
  xz.initialise_fftw_plan(n_threads, n_z, eh_vec);
  yx.initialise_fftw_plan(n_threads, n_x, eh_vec);
  yz.initialise_fftw_plan(n_threads, n_z, eh_vec);
  zx.initialise_fftw_plan(n_threads, n_x, eh_vec);
  zy.initialise_fftw_plan(n_threads, n_y, eh_vec);
}
