#include "hsl.h"

#define F77NAME(x) x ## _
// By default, the GNU Fortran compiler appends an underscore to the function names.
// If your Fortran compiler does not do this, you can change the definition of F77NAME accordingly.
// For example, one can set the following for Intel Fortran compiler:
//  "CMAKE_Fortran_FLAGS": "/names:lowercase /assume:underscore",
// to enable traditional Fortran name mangling.

extern "C"
{
    void F77NAME(mc29ad)(int* m, int* n, int* ne, double* a, int* ir, int* ic, double* r, double* c, double* w, int* lp, int* ifail);
    void F77NAME(mc75ad)(int* n, int* nz, int* la, double* a, int* ir, int* ic, double* cond, int* liw, int* iw, int* lw, double* w, int* icntl, int* info);
}

namespace Kratos
{
    void mc29ad_wrapper(int* m, int* n, int* ne, double* a, int* ir, int* ic, double* r, double* c, double* w, int* lp, int* ifail)
    {
        F77NAME(mc29ad)(m, n, ne, a, ir, ic, r, c, w, lp, ifail);
    }

    void mc75ad_wrapper(int* n, int* nz, int* la, double* a, int* ir, int* ic, double* cond, int* liw, int* iw, int* lw, double* w, int* icntl, int* info)
    {
        F77NAME(mc75ad)(n, nz, la, a, ir, ic, cond, liw, iw, lw, w, icntl, info);
    }
}
