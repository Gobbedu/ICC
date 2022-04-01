#ifndef __NEWTON_PADRAO__
#define __NEWTON_PADRAO__


#ifndef ____SIST_N_LINEAR__
#include "SistNlinear.h"
#endif

void NewtonPadrao(SistNl_t *snl, double *resposta, Tempo_t *time, int *nIter);

#endif
