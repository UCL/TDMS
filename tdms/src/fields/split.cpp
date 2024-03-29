#include "field.h"

#include <cstdlib>

#include "interpolation_methods.h"

using namespace std;


void SplitFieldComponent::initialise_fftw_plan(int n_threads, int size,
                                               EHVec &eh_vec) {

  this->n_threads = n_threads;
  plan_f = (fftw_plan *) malloc(sizeof(fftw_plan *) * n_threads);
  plan_b = (fftw_plan *) malloc(sizeof(fftw_plan *) * n_threads);

  for (int i = 0; i < n_threads; i++) {
    plan_f[i] = fftw_plan_dft_1d(size, eh_vec[i], eh_vec[i], FFTW_FORWARD,
                                 FFTW_MEASURE);
    plan_b[i] = fftw_plan_dft_1d(size, eh_vec[i], eh_vec[i], FFTW_BACKWARD,
                                 FFTW_MEASURE);
  }
}

void SplitFieldComponent::initialise_from_matlab(double ***tensor,
                                                 Dimensions &dims) {
  initialise(tensor, dims[2], dims[1], dims[0], true);
}

SplitFieldComponent::~SplitFieldComponent() {

  for (auto plan : {plan_f, plan_b}) {
    if (plan == nullptr) continue;

    for (int i = 0; i < n_threads; i++) { fftw_destroy_plan(plan[i]); }
    free(plan);
  }
}

void SplitField::allocate() {

  for (auto component : {&xy, &xz, &yx, &yz, &zx, &zy}) {
    component->allocate(tot.k + 1, tot.j + 1, tot.i + 1);
  }
}

SplitField::SplitField(int I_total, int J_total, int K_total) {
  tot.i = I_total;
  tot.j = J_total;
  tot.k = K_total;
}

void SplitField::zero() {

  for (auto component : {&xy, &xz, &yx, &yz, &zx, &zy}) { component->zero(); }
}

void SplitField::initialise_fftw_plan(int n_threads, EHVec &eh_vec) {

  // N_e_y
  int n_x = tot.i + delta_n() + 1;
  int n_y = tot.j + delta_n() + 1;
  int n_z = tot.k + delta_n() + 1;

  xy.initialise_fftw_plan(n_threads, n_y, eh_vec);
  xz.initialise_fftw_plan(n_threads, n_z, eh_vec);
  yx.initialise_fftw_plan(n_threads, n_x, eh_vec);
  yz.initialise_fftw_plan(n_threads, n_z, eh_vec);
  zx.initialise_fftw_plan(n_threads, n_x, eh_vec);
  zy.initialise_fftw_plan(n_threads, n_y, eh_vec);
}

double SplitField::largest_field_value() {
  double largest_value = 0.;
  for (int k = 0; k < (tot.k + 1); k++) {
    for (int j = 0; j < (tot.j + 1); j++) {
      for (int i = 0; i < (tot.i + 1); i++) {
        double x_field = fabs(xy(i, j, k) + xz(i, j, k));
        double y_field = fabs(yx(i, j, k) + yz(i, j, k));
        double z_field = fabs(zx(i, j, k) + zy(i, j, k));
        if (largest_value < x_field) { largest_value = x_field; }
        if (largest_value < y_field) { largest_value = y_field; }
        if (largest_value < z_field) { largest_value = z_field; }
      }
    }
  }
  return largest_value;
}
