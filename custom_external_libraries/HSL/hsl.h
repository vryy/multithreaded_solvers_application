#ifndef _HSL_H_
#define _HSL_H_

namespace Kratos
{
    void mc29ad_wrapper(int* m, int* n, int* ne, double* a, int* ir, int* ic, double* r, double* c, double* w, int* lp, int* ifail);
    void mc75ad_wrapper(int* n, int* nz, int* la, double* a, int* ir, int* ic, double* cond, int* liw, int* iw, int* lw, double* w, int* icntl, int* info);
}

#endif

