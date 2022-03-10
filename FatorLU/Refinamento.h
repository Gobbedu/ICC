#ifndef __REFINAMENTO__
#define __REFINAMENTO__

#include "SistLinear.h"

// aplica MAXIT iterações do Refinamento em SL com solução X
// EDITADO NA MACALAN
void refinamento(SistLinear_t *SL, double *X, int MAXIT);

void refinamentoLU(SistLinear_t *SL, double *X, int MAXIT);

void FatorLU(SistLinear_t *A, SistLinear_t *LU);

#endif
