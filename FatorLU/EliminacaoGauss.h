#ifndef __ELIM_GAUSS__
#define __ELIM_GAUSS__

#include "SistLinear.h"

void pivot(SistLinear_t *SL, int i);

void retrossubs(SistLinear_t *SL, double *X);

void triang(SistLinear_t *SL);

void eliminacaoGauss(SistLinear_t *SL, double *X);

#endif
