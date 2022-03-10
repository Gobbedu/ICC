#ifndef __REFINAMENTO__
#define __REFINAMENTO__

#include "SistLinear.h"

// aplica MAXIT iterações do Refinamento em SL com solução X
void refinamento(SistLinear_t *SL, double *X, int MAXIT);

#endif
