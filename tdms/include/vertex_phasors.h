/**
 * @file vertex_phasors.h
 * @brief Contains a class that handles the complex amplitude extraction at the
 * vertices.
 */
#pragma once

#include <complex>

#include "arrays/tdms_matrix.h"
#include "arrays/vector_typedefs.h"
#include "field.h"
#include "grid_labels.h"
#include "utils.h"

/**
 * Class container for handling complex amplitude samples at the vertices, and
 * their extraction.
 *
 * Abbreviated to CAmpSample in MATLAB code
 */
class VertexPhasors {
private:
  /* n_vertices()-by-3 int array.
  Each "row" corresponds to the index of a vertex at which to extract the field
  components requested.
  */
  Vertices vertices;
  /*! The MATLAB indices of the field components we want to extract at the
   * vertices.
   */
  FieldComponentsVector components;

  mxArray *mx_camplitudes = nullptr;//< Complex amplitudes at the vertices
  int f_ex_vector_size =
          0;//< Number of elements in the frequency extraction vector

  /* Storage for real and imag parts of mx_surface_amplitudes (these can be
  f_ex_vector_size * n_surface_vertices arrays of FullFieldSnapshots when MATLAB
  is removed!)

  Arrays are index by
  [frequency_index][field_component_position][vertex_id/number].

  frequency_index corresponds to the frequencies at which the user has requested
  we extract the amplitudes.

  field_component_position corresponds to components.index(field_component) of
  the field_component we are interested it. Since the user might not request
  consecutive components for extraction (EG [Ex, Ey, Hz] would correspond to [1,
  2, 6]) we need to be able to convert these MATLAB indices -> consecutive
  indices for storage here.
  */
  double ***camplitudesR = nullptr,
         ***camplitudesI = nullptr;//!< @copydoc camplitudesR

public:
  VertexPhasors() = default;
  /*! @copydoc set_from */
  VertexPhasors(const mxArray *ptr) { set_from(ptr); }

  /**
   * @brief Setup using data from an input file
   *
   * @param ptr Pointer to the struct containing the list of vertices and
   * components to extract phasors at/for
   */
  void set_from(const mxArray *ptr);

  /** @brief Get the pointer to the data */
  mxArray *get_mx_camplitudes() { return mx_camplitudes; }

  /**
   * @brief Allocate memory for the camplitude{R,I} arrays.
   *
   * Provided there are vertices for us to extract at, allocates the memory for
   * the camplitude{R,I} arrays and creates the mx_camplitudes pointer to the
   * output data.
   *
   * @param n_frequencies The number of frequencies at which we need to extract
   * phasors
   */
  void setup_complex_amplitude_arrays(int n_frequencies);

  /** @brief Fetch the number of vertices at which we are extracting phasors */
  int n_vertices() { return vertices.n_vertices(); }
  /** @brief Fetch the number of field components we are extracting */
  int n_components() { return components.size(); }
  /** @brief Returns true/false based on whether there are/aren't vertices to
   * extract at */
  bool there_are_vertices_to_extract_at() { return (n_vertices() > 0); }
  // Returns true/false based on whether there are/aren't elements in BOTH the
  // vertices and components arrays
  bool there_are_elements_in_arrays() {
    return (vertices.has_elements() &&
            tdms_vector_utils::has_elements(components));
  }

  /**
   * @brief Normalise the surface amplitudes at frequency_vector_index by the E-
   * and H-norms provided.
   *
   * E-field components in camplitudes{R,I} are divided by the (complex) Enorm.
   * H-field components in camplitudes{R,I} are divided by the (complex) Hnorm.
   *
   * @param frequency_index Frequency index, camplitudes{R,I}[frequency_index]
   * will be normalised
   * @param Enorm,Hnorm The {E,H}-norm to normalise the {E,H}-components by
   */
  void normalise_vertices(int frequency_index, std::complex<double> Enorm,
                          std::complex<double> Hnorm);

  /**
   * @brief Extract the phasor values at the vertices on the surface, for the
   * given frequency index
   *
   * @param frequency_index The entries in camplitudes{R,I}[frequency_index]
   * will be written to
   * @param E,H The electric,magnetic split field
   * @param n Current timestep index
   * @param omega Angular frequency
   * @param params The parameters for this simulation
   */
  void extractPhasorsVertices(int frequency_index, ElectricSplitField &E,
                              MagneticSplitField &H, int n, double omega,
                              SimulationParameters &params);

  /**
   * @brief Incriments camplitudes{R,I} at the given index by the field values
   * provided.
   *
   * If we allow element-wise assignment, we are essentially performing the
   * operations: camplitudesR[frequency_index][:][vertex_index] += F.real(),
   * camplitudesI[frequency_index][:][vertex_index] += F.imag().
   *
   * Only field components that we are extracting at the vertices are pulled out
   * of F.
   *
   * @param frequency_index Frequency index
   * @param vertex_index Vertex index
   * @param F Field values to assign
   */
  void update_vertex_camplitudes(int frequency_index, int vertex_index,
                                 FullFieldSnapshot F);

  ~VertexPhasors();
};
