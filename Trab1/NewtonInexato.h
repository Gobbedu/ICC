#include "utils.h"
#include "SistNlinear.h"
#include "NewtonPadrao.h"

void gauss_seidel(SistLinear_t *SL, double *X);
void NewtonInexato(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter);
