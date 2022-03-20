#ifndef __REFINAMENTO__
#define __REFINAMENTO__

#include "SistLinear.h"

// aplica MAXIT iterações do Refinamento em SL com solução X

//////////////////// FATORACAO LU ////////////////////
void refinamento(SistLinear_t *SL, double *X, int MAXIT);
void refinamentoLU(SistLinear_t *SL, double *X, int MAXIT);
void normsubs(SistLinear_t *SL, double *X);
void FatorLU(SistLinear_t *A, SistLinear_t *LU);

#endif
