#include "iterator_main_loop.h"

void IteratorMainLoop::log_update(ElectricSplitField &E_split, MagneticSplitField &H_split, int tind) {
  // if it has been long enough since we last printed the current iteration number
  if ((((double) time(NULL)) - last_logged_iteration_time) > seconds_between_logs) {

    double maxfield = std::max(E_split.largest_field_value(), H_split.largest_field_value());

    spdlog::info("Iterating: tind = {0:d}, maxfield = {1:e}", tind, maxfield);
    last_logged_iteration_time = double(time(NULL));
  }
}
