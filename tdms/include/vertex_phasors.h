/**
 * @file vertex_phasors.h
 * @brief Contains a class that handles the complex amplitude extraction at the vertices.
 */
#pragma once

#include <complex>

#include "arrays.h"
#include "field.h"
#include "grid_labels.h"

//Check if ComplexAmplitudeSample in array.h is redundant or not after making this class
class VertexPhasors {
private:
    Vertices vertices;
    FieldComponentsVector components;

    mxArray *mx_camplitudes = nullptr;//< Complex amplitudes at the vertices
    int f_ex_vector_size = 0;//< Number of elements in the frequency extraction vector

    /* Storage for real and imag parts of mx_surface_amplitudes (these can be f_ex_vector_size * n_surface_vertices arrays of FullFieldSnapshots when MATLAB is removed!)

    Arrays are index by [frequency_index][field component][vertex_id/number]. Frequency index corresponds to the frequencies at which the user has requested we extract the amplitudes.
    */
    double ***camplitudesR = nullptr, ***camplitudesI = nullptr;

public:
    VertexPhasors(const mxArray *ptr);

    void setup_camplitude_arrays(int _f_ex_vector_size);

    mxArray *get_mx_camplitudes() { return mx_camplitudes; }

    void normalise_vertices(int frequency_index, std::complex<double> Enorm, std::complex<double> Hnorm);

    int n_vertices(){ return vertices.n_vertices(); }
    bool there_are_vertices_to_extract_at() { return (n_vertices() > 0); }

    void extractPhasorsVertices(int frequency_index, ElectricSplitField &E, MagneticSplitField &H,
                                int n, double omega, SimulationParameters &params);

    void update_vertex_camplitudes(int frequency_index, int vertex_index, FullFieldSnapshot F);

    ~VertexPhasors();
};
