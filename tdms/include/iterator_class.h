/**
 * @file iterator_class.h
 * @brief Class that runs the TDMS simulation from being passed the inputs and command line arguments, running the FDTD solver, and passing the outputs back to the main() function for writing
 */
#pragma once

#include <spdlog/spdlog.h>

#include "matrix.h"

#include "iterator_executor.h"

/**
 * @brief Class that handles the setup, execution, output passing, and tear-down of the FDTD/PSTD simulation.
 *
 * This class is the lowest in the Iterator_ hierarchy, and explicitly defines methods for passing the results of the simulation back to the main() function via the plhs[] pointer-array.
 * It also handles the tear-down/cleanup of any malloc'd memory THAT IS NOT REQUIRED by the outputs themselves, but must be reserved when writing the outputs themselves (MATLAB is weird, TODO: Update when HDF5 comes in).
 *
 * This class handles the functionality of what was previously the file iterator.cpp.
 */
class Iterator : public Iterator_Executor {
public:
  Iterator(InputMatrices matrices_from_input_file, SolverMethod _solver_method,
           PreferredInterpolationMethods interpolation_method)
      : Iterator_Executor(matrices_from_input_file, _solver_method, interpolation_method){};

  /**
   * @brief Normalise the {E,H}-field phasors in the volume (if extracting phasors in the whole volume)
   */
  void normalise_field_volumes();
  /**
   * @brief Normalise the phasors on the user-defined surface (if extracting there)
   */
  void normalise_surface_phasors();
  /**
   * @brief Normalise the phasors at the user requested vertices (if there are any)
   */
  void normalise_vertex_phasors();
  /**
   * @brief Normalise the Id-array entries (if this is necessary)
   */
  void normalise_Id_arrays();

  void initialise_output_labels_from_input_labels() {
    output_grid_labels.initialise_from(input_grid_labels, E.il, E.iu, E.jl, E.ju, E.kl, E.ku);
  }

  /**
   * @brief Interpolate the extracted field values to the centres of the Yee cells, placing the interpolated values into the output.
   *
   * @param output_matrices The collection of output pointers
   */
  void interpolate_field_values(mxArray *output_matrices[]);
  /**
   * @brief Write the gridlabels (co-ordinates) of the interpolated fields to the output.
   *
   * @param output_matrices The collection of output pointers
   */
  void write_interpolated_gridlabels(mxArray *output_matrices[]);

  /**
   * @brief Regenerate the mesh of the user-defined surface.
   *
   * This is deleted prior to iterating to save memory, we recompute it after the main loop only if the user has requested it be exported.
   * Additionally, when we regenerate it we no longer want to re-initialise the surface_phasors we have computed, so we cannot recycle the setup_surface_mesh() method.
   *
   * @param output_matrices The collection of output pointers
   */
  void regenerate_mesh_for_facets(mxArray *output_matrices[]);

  /**
   * @brief Return the largest split-field value across the E- and H-split-fields
   */
  double compute_max_split_field_value() {
    return std::max(E_s.largest_field_value(), H_s.largest_field_value());
  }
};
