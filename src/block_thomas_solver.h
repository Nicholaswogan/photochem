#ifndef PHOTOCHEM_BLOCK_THOMAS_SOLVER_H
#define PHOTOCHEM_BLOCK_THOMAS_SOLVER_H

#include <sundials/sundials_linearsolver.h>

#ifdef __cplusplus
extern "C" {
#endif

SUNLinearSolver SUNLinSol_BlockThomas(N_Vector y, SUNMatrix A,
                                      sunindextype block_size);

#ifdef __cplusplus
}
#endif

#endif
