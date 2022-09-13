/**
 * @file fdtd_grid_initialiser.h
 * @brief Initialisation of the FDTD grid
 */
#pragma once
#include "mat_io.h"
#include <vector>
#include <string>

/**
 * @brief A class to initialise an FDTD grid from a MATLAB file
 */
class fdtdGridInitialiser {

private:
    const mxArray *pointer;         //< Pointer to the array
    const char *mat_filename;       //< Filename of the MATLAB file
    std::vector<mwSize> dimensions; //< The dimensions of the array

    /**
     * @brief Get a value from a integer attribute of the FDTD grid defined in a
     * .mat file
     * 
     * @param key the name of the attribute
     * @return int the value of the attribute
     */
    int value_of_attribute(const std::string& key);

public:
    /**
     * @brief Construct a new fdtd Grid Initialiser object
     * 
     * @param fdtd_pointer pointer to the FDTD grid
     * @param mat_filename the filename of the MATLAB file
     */
    fdtdGridInitialiser(const mxArray *fdtd_pointer, const char* mat_filename);

    /**
     * @brief Set an FDTD grid attribute to a tensor full of zeros
     * @param name of the attribute
     */
    void add_tensor(const std::string &name);
};
