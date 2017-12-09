#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP TDA_AlphaComplexDiag(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_AlphaComplexFiltration(SEXP, SEXP);
extern SEXP TDA_AlphaShapeDiag(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_AlphaShapeFiltration(SEXP, SEXP);
extern SEXP TDA_Bottleneck(SEXP, SEXP);
extern SEXP TDA_Dtm(SEXP, SEXP, SEXP);
extern SEXP TDA_DtmWeight(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_FiltrationDiag(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_FunFiltration(SEXP, SEXP);
extern SEXP TDA_GridDiag(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_GridFiltration(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_Kde(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_KdeDist(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_RipsDiag(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_RipsFiltration(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TDA_Wasserstein(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  { "TDA_AlphaComplexDiag",       (DL_FUNC)&TDA_AlphaComplexDiag,       5 },
  { "TDA_AlphaComplexFiltration", (DL_FUNC)&TDA_AlphaComplexFiltration, 2 },
  { "TDA_AlphaShapeDiag",         (DL_FUNC)&TDA_AlphaShapeDiag,         5 },
  { "TDA_AlphaShapeFiltration",   (DL_FUNC)&TDA_AlphaShapeFiltration,   2 },
  { "TDA_Bottleneck",             (DL_FUNC)&TDA_Bottleneck,             2 },
  { "TDA_Dtm",                    (DL_FUNC)&TDA_Dtm,                    3 },
  { "TDA_DtmWeight",              (DL_FUNC)&TDA_DtmWeight,              5 },
  { "TDA_FiltrationDiag",         (DL_FUNC)&TDA_FiltrationDiag,         5 },
  { "TDA_FunFiltration",          (DL_FUNC)&TDA_FunFiltration,          2 },
  { "TDA_GridDiag",               (DL_FUNC)&TDA_GridDiag,               7 },
  { "TDA_GridFiltration",         (DL_FUNC)&TDA_GridFiltration,         5 },
  { "TDA_Kde",                    (DL_FUNC)&TDA_Kde,                    6 },
  { "TDA_KdeDist",                (DL_FUNC)&TDA_KdeDist,                5 },
  { "TDA_RipsDiag",               (DL_FUNC)&TDA_RipsDiag,               8 },
  { "TDA_RipsFiltration",         (DL_FUNC)&TDA_RipsFiltration,         6 },
  { "TDA_Wasserstein",            (DL_FUNC)&TDA_Wasserstein,            3 },
  { NULL, NULL, 0 }
};

void R_init_TDA(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
