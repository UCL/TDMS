#include "timer.h"

#include <cstdio>

#include <omp.h>
#include <spdlog/spdlog.h>

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
  spdlog::info("Time ellapsed (s): %.03e", delta_seconds());
  start();
}
