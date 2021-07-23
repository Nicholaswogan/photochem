#include <stdio.h>
#include "cminpack.h"

void hybrd1c_wrapper(int (*fcn)(void *p, int n, const double *x, double *fvec, int iflag), 
						void* p, int* n, double* x, double* fvec, double* tol, int* info, double* wa,
						int* lwa)
{
  *info = __cminpack_func__(hybrd1)(fcn, p, *n, x, fvec, *tol, wa, *lwa);
}