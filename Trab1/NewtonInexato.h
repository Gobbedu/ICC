/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./testaSNL < sistemas.dat
********************************************************/

#include "utils.h"
#include "SistNlinear.h"
#include "NewtonPadrao.h"

void gauss_seidel(SistLinear_t *SL, double *X);
void NewtonInexato(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter);
