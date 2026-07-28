#include <iostream>
#include <string>
#include <cmath>
#include "validation.hpp"

int main()
{
    ::npop = 1; // validation expects global npop to be meaningful only for message text; pass 1
    int np = 1;
    int nsnp = 3;

    // allocate ld_in
    LDTYPE*** ld_in = new LDTYPE**[np];
    ld_in[0] = new LDTYPE*[nsnp];
    for(int j = 0; j < nsnp; ++j)
    {
        ld_in[0][j] = new LDTYPE[nsnp];
        for(int k = 0; k < nsnp; ++k)
            ld_in[0][j][k] = (j == k) ? 1.0f : 0.05f;
    }

    std::string err;
    // valid case
    if(!validate_ld_matrices(np, nsnp, ld_in, err))
    {
        std::cerr << "Valid case: unexpected failure: " << err << std::endl;
        return 1;
    }

    // invalid diagonal
    ld_in[0][1][1] = 0.9f;
    if(validate_ld_matrices(np, nsnp, ld_in, err))
    {
        std::cerr << "Invalid-diag case: expected failure but got success" << std::endl;
        return 1;
    }
    if(err.find("invalid diagonal") == std::string::npos)
    {
        std::cerr << "Invalid-diag case: unexpected error message: " << err << std::endl;
        return 1;
    }

    // restore diag and introduce asymmetry
    ld_in[0][1][1] = 1.0f;
    ld_in[0][0][1] = 0.1f;
    ld_in[0][1][0] = 0.2f;
    if(validate_ld_matrices(np, nsnp, ld_in, err))
    {
        std::cerr << "Asymmetry case: expected failure but got success" << std::endl;
        return 1;
    }
    if(err.find("not symmetric") == std::string::npos)
    {
        std::cerr << "Asymmetry case: unexpected error message: " << err << std::endl;
        return 1;
    }

    // NaN case
    ld_in[0][0][1] = ld_in[0][1][0] = 0.05f;
    ld_in[0][2][2] = NAN;
    if(validate_ld_matrices(np, nsnp, ld_in, err))
    {
        std::cerr << "NaN case: expected failure but got success" << std::endl;
        return 1;
    }
    if(err.empty())
    {
        std::cerr << "NaN case: expected an error message but got empty" << std::endl;
        return 1;
    }

    // cleanup
    for(int j = 0; j < nsnp; ++j)
        delete [] ld_in[0][j];
    delete [] ld_in[0];
    delete [] ld_in;

    std::cout << "validate_ld_matrices unit test: OK" << std::endl;
    return 0;
}
