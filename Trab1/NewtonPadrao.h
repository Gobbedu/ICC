#ifndef __NEWTON_PADRAO__
#define __NEWTON_PADRAO__

#include "utils.h"
#include "SistLinear.h"
#include "SistNlinear.h"

void NewtonPadrao(SistNl_t *snl, double *resposta, Tempo_t *time, int *nIter);

//////////////////// ELIMINACAO DE GAUSS ////////////////////
void pivot(SistLinear_t *SL, int i);
void triang(SistLinear_t *SL);
void retrossubs(SistLinear_t *SL, double *X);
void eliminacaoGauss(SistLinear_t *SL, double *X);


#endif
