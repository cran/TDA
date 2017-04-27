#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP TDA_AlphaComplexDiagGUDHI(SEXP, SEXP);
extern SEXP TDA_AlphaShapeDiagGUDHI(SEXP, SEXP);
extern SEXP TDA_Bottleneck(SEXP, SEXP);
extern SEXP TDA_Dtm(SEXP, SEXP, SEXP);
extern SEXP TDA_DtmWeight(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_GridDiag(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_Kde(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_KdeDist(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_RipsDiag(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_Wasserstein(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"TDA_AlphaComplexDiagGUDHI", (DL_FUNC) &TDA_AlphaComplexDiagGUDHI, 2},
    {"TDA_AlphaShapeDiagGUDHI",   (DL_FUNC) &TDA_AlphaShapeDiagGUDHI,   2},
    {"TDA_Bottleneck",            (DL_FUNC) &TDA_Bottleneck,            2},
    {"TDA_Dtm",                   (DL_FUNC) &TDA_Dtm,                   3},
    {"TDA_DtmWeight",             (DL_FUNC) &TDA_DtmWeight,             5},
    {"TDA_GridDiag",              (DL_FUNC) &TDA_GridDiag,              7},
    {"TDA_Kde",                   (DL_FUNC) &TDA_Kde,                   5},
    {"TDA_KdeDist",               (DL_FUNC) &TDA_KdeDist,               5},
    {"TDA_RipsDiag",              (DL_FUNC) &TDA_RipsDiag,              7},
    {"TDA_Wasserstein",           (DL_FUNC) &TDA_Wasserstein,           3},
    {NULL, NULL, 0}
};

void R_init_TDA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
