/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "utils.h"
#include "SistNlinear.h"
#include "NewtonPadrao.h"


void gauss_seidel(SistLinear_t *SL, double *X);
void NewtonInexato(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter);
