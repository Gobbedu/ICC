#ifndef __FATORACAO_LU__
#define __FATORACAO_LU__

#ifndef __SIST_LINEAR__
#include "SistLinear.h"
#endif

//////////////////// FATORACAO LU ////////////////////

void FatorLU(SistLinear_t *LU);
void EliminacaoLU(SistLinear_t *LU, double *X);
void normsubs(SistLinear_t *SL, double *X);

#endif
