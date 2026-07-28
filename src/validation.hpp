#pragma once
#include <string>
#include "data.hpp"

// Validate LD matrices: returns true if valid, false and sets err_msg otherwise.
bool validate_ld_matrices(int npop, int nsnp, LDTYPE*** ld_in, std::string &err_msg, float tol = 1e-6f);
