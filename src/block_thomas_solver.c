#include "block_thomas_solver.h"

#include <stdlib.h>
#include <string.h>

#include <sunmatrix/sunmatrix_band.h>
#include <sundials/sundials_dense.h>

#define ZERO RCONST(0.0)

typedef struct {
  sunindextype n;
  sunindextype block_size;
  sunindextype nblocks;
  sunindextype last_flag;
  realtype* diagonal;
  realtype* upper;
  realtype* lower;
  realtype* matrix_work;
  realtype* rhs_work;
  realtype* vector_work;
  realtype** columns;
  sunindextype* pivots;
} *SUNLinearSolverContent_BlockThomas;

#define BLOCK_CONTENT(S) ((SUNLinearSolverContent_BlockThomas)((S)->content))

static SUNLinearSolver_Type SUNLinSolGetType_BlockThomas(SUNLinearSolver S);
static SUNLinearSolver_ID SUNLinSolGetID_BlockThomas(SUNLinearSolver S);
static int SUNLinSolInitialize_BlockThomas(SUNLinearSolver S);
static int SUNLinSolSetup_BlockThomas(SUNLinearSolver S, SUNMatrix A);
static int SUNLinSolSolve_BlockThomas(SUNLinearSolver S, SUNMatrix A,
                                     N_Vector x, N_Vector b, realtype tol);
static sunindextype SUNLinSolLastFlag_BlockThomas(SUNLinearSolver S);
static int SUNLinSolSpace_BlockThomas(SUNLinearSolver S, long int* lenrwLS,
                                     long int* leniwLS);
static int SUNLinSolFree_BlockThomas(SUNLinearSolver S);

static void set_block_columns(SUNLinearSolverContent_BlockThomas content,
                              sunindextype block)
{
  sunindextype column;
  realtype* data = content->diagonal +
                   block*content->block_size*content->block_size;
  for (column = 0; column < content->block_size; column++) {
    content->columns[column] = data + column*content->block_size;
  }
}

static int matrix_has_block_tridiagonal_structure(
    SUNLinearSolverContent_BlockThomas content, SUNMatrix A)
{
  sunindextype column, row, row_min, row_max;
  sunindextype block_size = content->block_size;
  sunindextype n = content->n;
  sunindextype lower_bandwidth = SUNBandMatrix_LowerBandwidth(A);
  sunindextype upper_bandwidth = SUNBandMatrix_UpperBandwidth(A);

  for (column = 0; column < n; column++) {
    row_min = column > upper_bandwidth ? column - upper_bandwidth : 0;
    row_max = column + lower_bandwidth < n ?
              column + lower_bandwidth : n - 1;
    for (row = row_min; row <= row_max; row++) {
      if (row/block_size == column/block_size) continue;
      if (row%block_size == column%block_size) continue;
      if (SM_ELEMENT_B(A, row, column) != ZERO) return 0;
    }
  }
  return 1;
}

SUNLinearSolver SUNLinSol_BlockThomas(N_Vector y, SUNMatrix A,
                                      sunindextype block_size)
{
  SUNLinearSolver S = NULL;
  SUNLinearSolverContent_BlockThomas content = NULL;
  sunindextype n, nblocks, block_entries, coupling_entries;

  if (y == NULL || A == NULL || block_size <= 0) return NULL;
  if (SUNMatGetID(A) != SUNMATRIX_BAND) return NULL;
  if (SUNBandMatrix_Rows(A) != SUNBandMatrix_Columns(A)) return NULL;

  n = SUNBandMatrix_Rows(A);
  if (n <= 0 || n != N_VGetLength(y) || n%block_size != 0) return NULL;
  if (SUNBandMatrix_UpperBandwidth(A) < block_size ||
      SUNBandMatrix_LowerBandwidth(A) < block_size) return NULL;

  nblocks = n/block_size;
  block_entries = nblocks*block_size*block_size;
  coupling_entries = nblocks > 1 ? (nblocks - 1)*block_size : 0;

  S = SUNLinSolNewEmpty();
  if (S == NULL) return NULL;

  S->ops->gettype = SUNLinSolGetType_BlockThomas;
  S->ops->getid = SUNLinSolGetID_BlockThomas;
  S->ops->initialize = SUNLinSolInitialize_BlockThomas;
  S->ops->setup = SUNLinSolSetup_BlockThomas;
  S->ops->solve = SUNLinSolSolve_BlockThomas;
  S->ops->lastflag = SUNLinSolLastFlag_BlockThomas;
  S->ops->space = SUNLinSolSpace_BlockThomas;
  S->ops->free = SUNLinSolFree_BlockThomas;

  content = calloc(1, sizeof *content);
  if (content == NULL) {
    SUNLinSolFree(S);
    return NULL;
  }
  S->content = content;

  content->n = n;
  content->block_size = block_size;
  content->nblocks = nblocks;
  content->diagonal = malloc(block_entries*sizeof(realtype));
  content->upper = malloc(coupling_entries*sizeof(realtype));
  content->lower = malloc(coupling_entries*sizeof(realtype));
  content->matrix_work = malloc(block_size*block_size*sizeof(realtype));
  content->rhs_work = malloc(n*sizeof(realtype));
  content->vector_work = malloc(block_size*sizeof(realtype));
  content->columns = malloc(block_size*sizeof(realtype*));
  content->pivots = malloc(n*sizeof(sunindextype));

  if (content->diagonal == NULL ||
      (coupling_entries > 0 &&
       (content->upper == NULL || content->lower == NULL)) ||
      content->matrix_work == NULL || content->rhs_work == NULL ||
      content->vector_work == NULL || content->columns == NULL ||
      content->pivots == NULL) {
    SUNLinSolFree(S);
    return NULL;
  }

  content->last_flag = SUNLS_SUCCESS;
  return S;
}

static SUNLinearSolver_Type SUNLinSolGetType_BlockThomas(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_DIRECT;
}

static SUNLinearSolver_ID SUNLinSolGetID_BlockThomas(SUNLinearSolver S)
{
  return SUNLINEARSOLVER_CUSTOM;
}

static int SUNLinSolInitialize_BlockThomas(SUNLinearSolver S)
{
  if (S == NULL || S->content == NULL) return SUNLS_MEM_NULL;
  BLOCK_CONTENT(S)->last_flag = SUNLS_SUCCESS;
  return SUNLS_SUCCESS;
}

static int SUNLinSolSetup_BlockThomas(SUNLinearSolver S, SUNMatrix A)
{
  SUNLinearSolverContent_BlockThomas content;
  sunindextype block, row, column, info;
  sunindextype q;
  realtype* current;
  realtype* previous_upper;
  realtype* previous_lower;

  if (S == NULL || S->content == NULL || A == NULL) return SUNLS_MEM_NULL;
  if (SUNMatGetID(A) != SUNMATRIX_BAND) return SUNLS_ILL_INPUT;

  content = BLOCK_CONTENT(S);
  q = content->block_size;
  if (SUNBandMatrix_Rows(A) != content->n ||
      SUNBandMatrix_Columns(A) != content->n ||
      SUNBandMatrix_UpperBandwidth(A) < q ||
      SUNBandMatrix_LowerBandwidth(A) < q ||
      !matrix_has_block_tridiagonal_structure(content, A)) {
    content->last_flag = SUNLS_ILL_INPUT;
    return SUNLS_ILL_INPUT;
  }

  for (block = 0; block < content->nblocks; block++) {
    current = content->diagonal + block*q*q;
    for (column = 0; column < q; column++) {
      for (row = 0; row < q; row++) {
        current[column*q + row] =
            SM_ELEMENT_B(A, block*q + row, block*q + column);
      }
    }
  }
  for (block = 0; block + 1 < content->nblocks; block++) {
    for (row = 0; row < q; row++) {
      content->upper[block*q + row] =
          SM_ELEMENT_B(A, block*q + row, (block + 1)*q + row);
      content->lower[block*q + row] =
          SM_ELEMENT_B(A, (block + 1)*q + row, block*q + row);
    }
  }

  set_block_columns(content, 0);
  info = denseGETRF(content->columns, q, q, content->pivots);
  if (info != 0) {
    content->last_flag = info;
    return SUNLS_LUFACT_FAIL;
  }

  for (block = 1; block < content->nblocks; block++) {
    previous_upper = content->upper + (block - 1)*q;
    previous_lower = content->lower + (block - 1)*q;
    memset(content->matrix_work, 0, q*q*sizeof(realtype));
    for (column = 0; column < q; column++) {
      content->matrix_work[column*q + column] = previous_upper[column];
      denseGETRS(content->columns, q, content->pivots + (block - 1)*q,
                 content->matrix_work + column*q);
    }

    current = content->diagonal + block*q*q;
    for (column = 0; column < q; column++) {
      for (row = 0; row < q; row++) {
        current[column*q + row] -=
            previous_lower[row]*content->matrix_work[column*q + row];
      }
    }

    set_block_columns(content, block);
    info = denseGETRF(content->columns, q, q,
                      content->pivots + block*q);
    if (info != 0) {
      content->last_flag = info;
      return SUNLS_LUFACT_FAIL;
    }
  }

  content->last_flag = SUNLS_SUCCESS;
  return SUNLS_SUCCESS;
}

static int SUNLinSolSolve_BlockThomas(SUNLinearSolver S, SUNMatrix A,
                                     N_Vector x, N_Vector b, realtype tol)
{
  SUNLinearSolverContent_BlockThomas content;
  sunindextype block, row, q;
  realtype* xdata;
  realtype* bdata;
  realtype* rhs;

  (void)A;
  (void)tol;
  if (S == NULL || S->content == NULL || x == NULL || b == NULL)
    return SUNLS_MEM_NULL;

  content = BLOCK_CONTENT(S);
  q = content->block_size;
  if (N_VGetLength(x) != content->n || N_VGetLength(b) != content->n) {
    content->last_flag = SUNLS_ILL_INPUT;
    return SUNLS_ILL_INPUT;
  }
  xdata = N_VGetArrayPointer(x);
  bdata = N_VGetArrayPointer(b);
  if (xdata == NULL || bdata == NULL) {
    content->last_flag = SUNLS_MEM_FAIL;
    return SUNLS_MEM_FAIL;
  }

  rhs = content->rhs_work;
  memcpy(rhs, bdata, content->n*sizeof(realtype));

  for (block = 1; block < content->nblocks; block++) {
    memcpy(content->vector_work, rhs + (block - 1)*q,
           q*sizeof(realtype));
    set_block_columns(content, block - 1);
    denseGETRS(content->columns, q, content->pivots + (block - 1)*q,
               content->vector_work);
    for (row = 0; row < q; row++) {
      rhs[block*q + row] -=
          content->lower[(block - 1)*q + row]*content->vector_work[row];
    }
  }

  block = content->nblocks - 1;
  memcpy(xdata + block*q, rhs + block*q, q*sizeof(realtype));
  set_block_columns(content, block);
  denseGETRS(content->columns, q, content->pivots + block*q,
             xdata + block*q);

  while (block > 0) {
    block--;
    for (row = 0; row < q; row++) {
      xdata[block*q + row] = rhs[block*q + row] -
          content->upper[block*q + row]*xdata[(block + 1)*q + row];
    }
    set_block_columns(content, block);
    denseGETRS(content->columns, q, content->pivots + block*q,
               xdata + block*q);
  }

  content->last_flag = SUNLS_SUCCESS;
  return SUNLS_SUCCESS;
}

static sunindextype SUNLinSolLastFlag_BlockThomas(SUNLinearSolver S)
{
  if (S == NULL || S->content == NULL) return -1;
  return BLOCK_CONTENT(S)->last_flag;
}

static int SUNLinSolSpace_BlockThomas(SUNLinearSolver S, long int* lenrwLS,
                                     long int* leniwLS)
{
  SUNLinearSolverContent_BlockThomas content;
  sunindextype q, coupling_entries;
  if (S == NULL || S->content == NULL) return SUNLS_MEM_NULL;
  if (lenrwLS == NULL || leniwLS == NULL) return SUNLS_ILL_INPUT;
  content = BLOCK_CONTENT(S);
  q = content->block_size;
  coupling_entries = content->nblocks > 1 ?
                     (content->nblocks - 1)*q : 0;
  *lenrwLS = (long int)(content->nblocks*q*q + 2*coupling_entries +
                        q*q + content->n + q);
  *leniwLS = (long int)content->n;
  return SUNLS_SUCCESS;
}

static int SUNLinSolFree_BlockThomas(SUNLinearSolver S)
{
  SUNLinearSolverContent_BlockThomas content;
  if (S == NULL) return SUNLS_SUCCESS;
  content = BLOCK_CONTENT(S);
  if (content != NULL) {
    free(content->diagonal);
    free(content->upper);
    free(content->lower);
    free(content->matrix_work);
    free(content->rhs_work);
    free(content->vector_work);
    free(content->columns);
    free(content->pivots);
    free(content);
    S->content = NULL;
  }
  if (S->ops != NULL) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  return SUNLS_SUCCESS;
}
