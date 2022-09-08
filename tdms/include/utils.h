/**
 * @file utils.h
 *
 * Defines utility functions.
 */

/**
 * @brief Throws a runtime error if a file is not found.
 * 
 * @param filename The name of the file to check.
 * @param mode The mode to try and open with.
 */
void assert_can_open_file(const char* filename, const char* mode);

/**
 * @brief Check two strings are equal
 * 
 * @param a The first string
 * @param b The second string
 * @return true if the strings are the same
 * @return false otherwise
 */
bool are_equal(const char* a, const char* b);
