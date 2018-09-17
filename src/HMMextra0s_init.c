#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(loop1)(int *m, int *T, double *phi, double *pRS, double *gamma, double *logalp, double *lscale, double *tmp);

extern void F77_NAME(loop2)(int *m, int *T, double *phi, double *pRS, double *gamma, double *logbet, double *lscale, double *tmp);

static const R_FortranMethodDef FortranEntries[] = {
    {"loop1", (DL_FUNC) &F77_NAME(loop1), 8},
    {"loop2", (DL_FUNC) &F77_NAME(loop2), 8},
    {NULL, NULL, 0}
};

void R_init_HMMextra0s(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
