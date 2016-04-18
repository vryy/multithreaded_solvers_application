#include "hsl.h"

namespace Kratos
{

    void mc29ad(int* m, int* n, int* ne, double* a, int* ir, int* ic, double* r, double* c, double* w, int* lp, int* ifail)
    {
        F77NAME(mc29ad)(m, n, ne, a, ir, ic, r, c, w, lp, ifail);
    }

    void mc75ad(int* n, int* nz, int* la, double* a, int* ir, int* ic, double* cond, int* liw, int* iw, int* lw, double* w, int* icntl, int* info)
    {
        F77NAME(mc75ad)(n, nz, la, a, ir, ic, cond, liw, iw, lw, w, icntl, info);
    }
}
