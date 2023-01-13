#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "utils.h"

/* .C calls */
extern void d011ch(
  double * z, 
  int * d, 
  int *dup, 
  double *sur, 
  double *jum, 
  int * max, 
  double * err, 
  int * r, 
  int * s, 
  int * rs, 
  double *zext, 
  int * dext, 
  double *wext,
  double * pK, 
  double * pkonst, 
  double *konstdist, 
  double *konstjump, 
  double * llratio);

extern void urnew010(
    double * z, int * d, int *dup, double * sur, double * jum, int * max, 
    double * err, int * r, int * s, int * rs, 
    double * zext, int * dext, double * wext);

static const R_CMethodDef CEntries[] = {
    {"d011ch",   (DL_FUNC) &d011ch,   18},
    {"urnew010", (DL_FUNC) &urnew010, 13},
    {NULL, NULL, 0}
};

void R_init_dblcens(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
