#ifndef _HSL_H_
#define _HSL_H_

#define F77NAME(x) x ## _

extern "C"
{
    void F77NAME(mc29ad)(int* m, int* n, int* ne, double* a, int* ir, int* ic, double* r, double* c, double* w, int* lp, int* ifail);
    void F77NAME(mc75ad)(int* n, int* nz, int* la, double* a, int* ir, int* ic, double* cond, int* liw, int* iw, int* lw, double* w, int* icntl, int* info);
}

namespace Kratos
{
    void mc29ad(int* m, int* n, int* ne, double* a, int* ir, int* ic, double* r, double* c, double* w, int* lp, int* ifail);
    void mc75ad(int* n, int* nz, int* la, double* a, int* ir, int* ic, double* cond, int* liw, int* iw, int* lw, double* w, int* icntl, int* info);
}

#endif

