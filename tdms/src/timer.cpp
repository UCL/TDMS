#include "cstdio"
#include "omp.h"
#include "timer.h"


void Timer::start() {
  start_time = omp_get_wtime();
}

void Timer::end() {
  end_time = omp_get_wtime();
}

double Timer::delta_seconds() const {
  return start_time - end_time;
}

/**
 * Log the difference in time and reset the timer
 */
void Timer::click(){
  end();
  fprintf(stdout, "∆time = %.03e s\n", delta_seconds());
  start();
}
