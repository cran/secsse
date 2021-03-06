#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(secsse_fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(secsse_initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(secsse_runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(secsse_runmod2)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(secsse_runmod_ct)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(cla_secsse_runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(cla_secsse_runmod_ct_e)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(cla_secsse_runmod_ct_d)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"secsse_fill1d", (DL_FUNC) &F77_NAME(secsse_fill1d), 4},
  {"secsse_initmod", (DL_FUNC) &F77_NAME(secsse_initmod), 1},
  {"secsse_runmod", (DL_FUNC) &F77_NAME(secsse_runmod), 6},
  {"secsse_runmod2", (DL_FUNC) &F77_NAME(secsse_runmod2), 6},
  {"secsse_runmod_ct", (DL_FUNC) &F77_NAME(secsse_runmod_ct), 6},
  {"cla_secsse_runmod", (DL_FUNC) &F77_NAME(cla_secsse_runmod), 6},
  {"cla_secsse_runmod_ct_e", (DL_FUNC) &F77_NAME(cla_secsse_runmod_ct_e), 6},
  {"cla_secsse_runmod_ct_d", (DL_FUNC) &F77_NAME(cla_secsse_runmod_ct_d), 6},
  {NULL, NULL, 0}
};

void R_init_secsse(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
