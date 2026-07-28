#include "validation.hpp"
#include <cmath>
#include <sstream>

bool validate_ld_matrices(int npop, int nsnp, LDTYPE*** ld_in, std::string &err_msg, float tol)
{
    if(npop <= 0)
    {
        err_msg = "npop must be > 0";
        return false;
    }
    for(int i = 0; i < npop; ++i)
    {
        for(int j = 0; j < nsnp; ++j)
        {
            float diag = ld_in[i][j][j];
            if(!std::isfinite(diag) || fabsf(diag - 1.0f) > tol)
            {
                std::ostringstream oss;
                oss << "LD matrix for population " << i << " has invalid diagonal at variant " << j << ": " << diag;
                err_msg = oss.str();
                return false;
            }
        }
        for(int j = 0; j < nsnp; ++j)
        {
            for(int k = j+1; k < nsnp; ++k)
            {
                float a = ld_in[i][j][k];
                float b = ld_in[i][k][j];
                if(!std::isfinite(a) || !std::isfinite(b) || fabsf(a - b) > tol)
                {
                    std::ostringstream oss;
                    oss << "LD matrix for population " << i << " is not symmetric at (" << j << "," << k << "): " << a << " vs " << b;
                    err_msg = oss.str();
                    return false;
                }
            }
        }
    }
    return true;
}
