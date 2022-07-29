#include "string"

/**
 * A collection of matlab matrices with names
 */
class MatrixCollection {

public:
    int n_matrices = 0;
    char** matrix_names = nullptr;

    MatrixCollection();
    explicit MatrixCollection(char **names, int number);

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
