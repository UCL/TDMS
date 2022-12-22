/**
 * @file matrix_collection.h
 * @brief A collection of named MATLAB matrices.
 */
#pragma once

#include <vector>
#include <string>

#include "mat_io.h"

/**
 * A collection of matlab matrices with names
 */
class MatrixCollection {

public:
    int n_matrices = 0;            //< The number of matrices in the collection
    std::vector<std::string> matrix_names; //< The names of the matrices

    /** @brief Construct a new Matrix Collection object */
    MatrixCollection() = default;

    /**
     * @brief Construct a new Matrix Collection object
     *
     * @param names of the matrices
     * @param number of matrices in the collection
     */
    explicit MatrixCollection(std::vector<std::string> names);

    /**
     * @brief Check we have as many or more matrices than another
     * MatrixCollection
     *
     * Throws a runtime error if there are fewer matrices in the other
     * MatrixCollection
     *
     * @param other the other MatrixCollection
     */
    void check_has_at_least_as_many_matrices_as(MatrixCollection &other);
};


/**
 * A collection of matlab matrices with names created from a .mat file
 */
class MatFileMatrixCollection : public MatrixCollection{

public:
    MATFile* mat_file;

    explicit MatFileMatrixCollection(const char *filename);

    ~MatFileMatrixCollection();
};
