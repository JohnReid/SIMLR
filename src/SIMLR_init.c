#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .Call calls */
extern SEXP projsplx(SEXP, SEXP);
extern SEXP _SIMLR_Rtsne_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"projsplx",         (DL_FUNC) &projsplx,          2},
  {"_SIMLR_Rtsne_cpp", (DL_FUNC) &_SIMLR_Rtsne_cpp, 10},
  {NULL, NULL, 0}
};

void R_init_SIMLR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

