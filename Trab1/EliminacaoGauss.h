#ifndef __ELIM_GAUSS__
#define __ELIM_GAUSS__

#include "utils.h"
#include "SistLinear.h"

//////////////////// ELIMINACAO DE GAUSS ////////////////////
void pivot(SistLinear_t *SL, int i);
void triang(SistLinear_t *SL);
void retrossubs(SistLinear_t *SL, double *X);
void eliminacaoGauss(SistLinear_t *SL, double *X);


#endif
