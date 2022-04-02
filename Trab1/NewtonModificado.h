#ifndef __FATORACAO_LU__
#define __FATORACAO_LU__

#include "utils.h"
#include "SistNlinear.h"
#include "EliminacaoGauss.h"

#define HESS_STEP 1 // usa # de variaveis

//////////////////// FATORACAO LU ////////////////////

void NewtonModificado(SistNl_t *snl, double *resposta, Tempo_t *tempo, int *numIteracoes);
void FatorLU(SistLinear_t *LU);
void EliminacaoLU(SistLinear_t *LU, double *X);
void normsubs(SistLinear_t *SL, double *X);

#endif
