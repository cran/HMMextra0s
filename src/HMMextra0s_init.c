#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(prsloop)(int *m, int *nn, double *pie, double *R, double *mu, double *sig, double *Z, double *pRS);

extern void F77_NAME(loop1)(int *m, int *T, double *phi, double *pRS, double *gamma, double *logalp, double *lscale, double *tmp);

extern void F77_NAME(loop2)(int *m, int *T, double *phi, double *pRS, double *gamma, double *logbet, double *lscale, double *tmp);

extern void F77_NAME(estep)(int *m, int *nn, double *logalpha, double *logbeta, double *ll, double *pRS, double *gamma, double *v, double *w);

extern void F77_NAME(mstep1d)(int *n, int *m, int *nn, double *v, double *Z, double *R, double *hatpie, double *hatmu, double *hatsig);

extern void F77_NAME(mstep2d)(int *n, int *m, int *nn, double *v, double *Z, double *R, double *hatpie, double *hatmu, double *hatsig);

static const R_FortranMethodDef FortranEntries[] = {
    {"prsloop", (DL_FUNC) &F77_NAME(prsloop), 8},
    {"loop1", (DL_FUNC) &F77_NAME(loop1), 8},
    {"loop2", (DL_FUNC) &F77_NAME(loop2), 8},
    {"estep", (DL_FUNC) &F77_NAME(estep), 9},
    {"mstep1d", (DL_FUNC) &F77_NAME(mstep1d), 9},
    {"mstep2d", (DL_FUNC) &F77_NAME(mstep2d), 9},
    {NULL, NULL, 0}
};

void R_init_HMMextra0s(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
