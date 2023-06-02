#include "vertex_phasors.h"

#include <spdlog/spdlog.h>

#include "globals.h"
#include "matlabio.h"
#include "matrix.h"

using namespace std;
using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;

void VertexPhasors::set_from(const mxArray *ptr) {
  // extract information from the MATLAB struct that was passed
  if (mxIsEmpty(ptr)) {
    spdlog::info("VertexPhasors: struct provided is empty");
    return;
  }
  assert_is_struct_with_n_fields(ptr, 2,
                                 "VertexPhasors (using campssample array)");
  vertices.initialise(ptr);

  // Setup the components by reading the MATLAB buffer
  // This should use the HDF5Reader methods when this class is detached from
  // MATLAB. This section will also need to use the to_vector_int method from
  // tdms_vector_utils.
  mxArray *components_field =
          ptr_to_matrix_in(ptr, "components", "campssample");
  if (!mxIsEmpty(components_field)) {
    const mwSize *dims = mxGetDimensions(components_field);
    unsigned int n_component_elements = max(dims[0], dims[1]);
    int *component_entries_as_doubles = (int *) mxGetPr(components_field);
    components.resize(n_component_elements);
    for (unsigned int c_index = 0; c_index < n_component_elements; c_index++) {
      components[c_index] = component_entries_as_doubles[c_index];
    }
  }
}

void VertexPhasors::normalise_vertices(int frequency_index,
                                       complex<double> Enorm,
                                       complex<double> Hnorm) {
  for (int i = FieldComponents::Ex; i <= FieldComponents::Hz; i++) {
    // determine whether we are extracting this field component at the vertices
    // or not
    int index_in_camplitudes = tdms_vector_utils::index(components, i);
    if (index_in_camplitudes >= 0) {
      // we are extracting this component, determine the normalisation factor
      complex<double> norm = i <= FieldComponents::Ez ? Enorm : Hnorm;
      // loop over all entries and normalise
      for (int vindex = 0; vindex < n_vertices(); vindex++) {
        complex<double> normalised_amplitude =
                complex<double>(camplitudesR[frequency_index]
                                            [index_in_camplitudes][vindex],
                                camplitudesI[frequency_index]
                                            [index_in_camplitudes][vindex]) /
                norm;
        camplitudesR[frequency_index][index_in_camplitudes][vindex] =
                normalised_amplitude.real();
        camplitudesI[frequency_index][index_in_camplitudes][vindex] =
                normalised_amplitude.imag();
      }
    }
  }
}

void VertexPhasors::setup_complex_amplitude_arrays(int n_frequencies) {
  // store the size of the frequency extraction vector
  f_ex_vector_size = n_frequencies;

  // setup camplitudes{R,I} storage if we need them for this run
  int dims[3] = {0, 0, 0};
  if (there_are_vertices_to_extract_at()) {
    dims[0] = n_vertices();
    dims[1] = components.size();
    dims[2] = f_ex_vector_size;
    mx_camplitudes = mxCreateNumericArray(3, (const mwSize *) dims,
                                          mxDOUBLE_CLASS, mxCOMPLEX);
    camplitudesR = cast_matlab_3D_array(mxGetPr(mx_camplitudes), dims[0],
                                        dims[1], dims[2]);
    camplitudesI = cast_matlab_3D_array(mxGetPi(mx_camplitudes), dims[0],
                                        dims[1], dims[2]);
  }
}

void VertexPhasors::extractPhasorsVertices(int frequency_index,
                                           ElectricSplitField &E,
                                           MagneticSplitField &H, int n,
                                           double omega,
                                           SimulationParameters &params) {
  int vindex;
  FullFieldSnapshot F;
  complex<double> phaseTermE, phaseTermH, cphaseTermE, cphaseTermH;

  phaseTermE = fmod(omega * ((double) n) * params.dt, 2 * DCPI);
  phaseTermH = fmod(omega * ((double) n + 0.5) * params.dt, 2 * DCPI);

  cphaseTermH = exp(phaseTermH * IMAGINARY_UNIT) * 1. / ((double) params.Npe);
  cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. / ((double) params.Npe);

  /* Loop over every vertex requested by the user.

  Since the value of the phasors at each vertex is entirely determined from the
  previously calculated fields, these computations can be done in parallel as
  the computation of the phasors is independent of one another. Ergo, we use a
  parallel loop.
  */
#pragma omp parallel default(shared) private(F, vindex)
  {
#pragma omp for
    for (vindex = 0; vindex < n_vertices(); vindex++) {// loop over every vertex
      CellCoordinate current_cell{vertices[0][vindex], vertices[1][vindex],
                                  vertices[2][vindex]};

      switch (params.dimension) {
        case Dimension::THREE:
          F.Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
          F.Ey = E.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          F.Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          F.Hx = H.interpolate_to_centre_of(AxialDirection::X, current_cell);
          F.Hy = H.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          F.Hz = H.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          break;
        case Dimension::TRANSVERSE_ELECTRIC:
          F.Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
          F.Ey = E.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          F.Hz = H.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          break;
        case Dimension::TRANSVERSE_MAGNETIC:
          F.Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
          F.Hx = H.interpolate_to_centre_of(AxialDirection::X, current_cell);
          F.Hy = H.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          break;
        default:
          throw runtime_error("Dimension was not recognised!");
          break;
      }
      // multiply by phasor factors
      F.multiply_E_by(cphaseTermE);
      F.multiply_H_by(cphaseTermH);

      // update the master arrays
      update_vertex_camplitudes(frequency_index, vindex, F);
    }
  }// end parallel region
}

void VertexPhasors::update_vertex_camplitudes(int frequency_index,
                                              int vertex_index,
                                              FullFieldSnapshot F) {
  for (int component_id = FieldComponents::Ex;
       component_id <= FieldComponents::Hz; component_id++) {
    int idx = tdms_vector_utils::index(components, component_id);
    // MATLAB indexes Ex->Hz with 1->6 (FieldComponents enum) whilst C++ indexes
    // Ex->Hz with 0->5 (FullFieldSnapshot) this is not an ideal fix, but it is
    // the simplest so long as we remember the correspondence above
    int cpp_component_index = component_id - 1;
    if (idx >= 0) {
      // this component has been requested for extraction
      camplitudesR[frequency_index][idx][vertex_index] +=
              real(F[cpp_component_index]);
      camplitudesI[frequency_index][idx][vertex_index] +=
              imag(F[cpp_component_index]);
    }
  }
}

VertexPhasors::~VertexPhasors() {
  // cleanup storage arrays if they were created
  if (there_are_vertices_to_extract_at()) {
    free_cast_matlab_3D_array(camplitudesR, f_ex_vector_size);
    free_cast_matlab_3D_array(camplitudesI, f_ex_vector_size);
  }
}
