/**
 * @file fieldsample.h
 * @brief Class that stores the extracted field values at user-defined vertices.
 */
#pragma once

#include "arrays.h"
#include "field.h"
#include "simulation_parameters.h"

class FieldSample{

private:
  double**** tensor = nullptr;

public:
  mxArray* mx;       // Matlab array

  Vector<int> i;     // Indices along the x-direction of locations at which to sample the field
  Vector<int> j;     // Indices along the y-direction of locations at which to sample the field
  Vector<int> k;     // Indices along the z-direction of locations at which to sample the field
  Vector<double> n;  // Vector of the moments of the field to sample

  explicit FieldSample(const mxArray *ptr);

  // Return true if all vectors in this instance are non-empty (have size > 0)
  bool all_vectors_are_non_empty() const{
          return i.size() > 0 && j.size() > 0 && k.size() > 0 && n.size() > 0;
  };

  inline double*** operator[] (int value) const { return tensor[value]; };

  /**
   * @brief Extract the (Electric) field values at the vertices
   *
   * @param E_split Values of the electric (split) field
   * @param pml A description of the perfectly matched layer being used in this simulation
   * @param n_simulation_timesteps The (total) number of timesteps in this simulation
   */
  void extract(ElectricSplitField E_split, PerfectlyMatchedLayer pml, int n_simulation_timesteps);

  ~FieldSample();
};
