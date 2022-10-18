#include "timer.h"

#include <cstdio>

#include <omp.h>


void Timer::start() {
  start_time = omp_get_wtime();
}

void Timer::end() {
  end_time = omp_get_wtime();
}

double Timer::delta_seconds() const {
  return end_time - start_time ;
}

void Timer::click(){
  end();
  fprintf(stdout, "∆time = %.03e s\n", delta_seconds());
  start();
}
