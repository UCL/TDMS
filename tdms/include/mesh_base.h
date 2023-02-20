/**
 * @file mesh_base.h
 * @brief Generation of orientated mesh.
 *
 * Generate an oriented mesh on the surface of a cuboid within the FDTD grid.
 */
#pragma once

#include "mat_io.h"
#include "simulation_parameters.h"

/**
 * @brief Generate a matrix of vertices which define a triangulation of a
 * regular two dimensional grid.
 *
 * This function assumes that the space of interest is a 2d surface with
 * coordinates (i,j). I0 represents the lowest value of i for any point on the
 * rectangular grid and I1 the highest. Similarly for j. A value of k is
 * constant. A line of the output matrix looks like:
 *
 *   i1 j1 k i2 j2 k i3 j3 k
 *
 * Triangles are taken by subdividing squares in the grid in a regular manner.
 *
 * @param coordmap is an integer array with three entries. This array can be a
 * permutation of {0,1,2}. This defines the mapping between i,j,k and the
 * indices in the output matrix.  For example, if coordmap = {0,1,2} then a row
 * in the matric would look like:
 *
 *   i1 j1 k i2 j2 k i3 j3 k
 *
 * If, however, we have coordmap = {2,1,0} then we would get
 *
 *   k j1 i1 k j2 i2 k j3 i3
 *
 * This should be interpreted as original i columns moves to column k. original
 * k column moves to column i.
 *
 *         i
 *      I1 ^ .  .  .  .
 *         | .  .  .  .
 *      I0 | .  .  .  .
 *         +------------>j
 *          J0        J1
 *
 * @param order specifies the direction of the surface normals of the triangles.
 * This can take only 2 possible values +1 or -1. They have the following
 * meaning:
 *
 * order = 1 means that the surface normal for a triangle in the:
 *     xy plane will || to the z-axis
 * 	   zy plane will || to the x-axis
 *     xz plane will || to the negative z-axis
 *
 * order = -1 means surface normals are in the opposite direction. The surface
 * normal is assumed to be in the direction (p2-p1)x(p3-p1) where p1-p3 are the
 * points which define the triangle, in the order that they are listed in the
 * facet matrix.
 *
 * The space allocated by *vertexMatrix must be freed after use.
 */
void triangulatePlane(int I0, int I1, int J0, int J1, int K, int coordmap[],
                      int order, mxArray **vertexMatrix);
void triangulatePlaneSkip(int I0, int I1, int J0, int J1, int K, int coordmap[],
                          int order, mxArray **vertexMatrix, int dI, int dJ);
/**
 * @brief
 *
 * @param vertexMatrix should be a 6 element array. Generates 6 arrays of facets
 * using triangulatePlane. Each matrix is a plane of the cuboid which is defined
 * by:
 *
 * (I0,I1)x(J0,J1)x(K0,K1)
 *
 * Each vertexMatrix[i] should be destroyed after calling this function
 */
void triangulateCuboid(int I0, int I1, int J0, int J1, int K0, int K1,
                       mxArray **vertexMatrix);
void triangulateCuboidSkip(int I0, int I1, int J0, int J1, int K0, int K1,
                           mxArray **vertexMatrix, int dI, int dJ, int dK);

/**
 * @brief Generates a triangulation of a cuboid defined the surface of a regular
 * grid.
 *
 * The result is returned in a concise manner, ie, a list of vertices and a list
 * of facets which index in to the list of vertices.
 *
 * The list of vertices is itself a list of indices in to the x, y and z grid
 * label vectors. In this sense this function deals only with the topology of
 * the cuboid and the mesh. An extra step is required to generate the actual
 * mesh from the values returned by this function.
 *
 * The surface of the volume [I0,I1]x[J0,J1]x[K0,K1] is meshed by this function.
 *
 * @param vertices an array of vertices, each row is a numbered vertex.
 * @param facets an array of facets each of which is created using 3 vertex
 * indices. Each index is an index in to the vertices array.
 */
void conciseTriangulateCuboid(int I0, int I1, int J0, int J1, int K0, int K1,
                              mxArray **vertices, mxArray **facets);
void conciseTriangulateCuboidSkip(int I0, int I1, int J0, int J1, int K0,
                                  int K1,
                                  const SurfaceSpacingStride &spacing_stride,
                                  mxArray **vertices, mxArray **facets);

void conciseCreateBoundary(int I0, int I1, int K0, int K1, mxArray **vertices,
                           mxArray **facets);
