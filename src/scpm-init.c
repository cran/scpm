#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(h2d)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"h2d", (DL_FUNC) &F77_NAME(h2d), 5},
    {NULL, NULL, 0}
};

void R_init_scpm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
