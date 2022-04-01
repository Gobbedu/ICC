#ifndef __NEWTON_PADRAO__
#define __NEWTON_PADRAO__


#include "SistNlinear.h"
#include "utils.h"
#include "EliminacaoGauss.h"

void NewtonPadrao(SistNl_t *snl, double *resposta, Tempo_t *time, int *nIter);

#endif
