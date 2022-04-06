/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./testaSNL < sistemas.dat
********************************************************/

#ifndef __NEWT_MODIFICADO__
#define __NEWT_MODIFICADO__

#include "utils.h"
#include "SistNlinear.h"
#include "NewtonPadrao.h"

#define HESS_STEP 1 // usa # de variaveis

//////////////////// FATORACAO LU ////////////////////

void NewtonModificado(SistNl_t *snl, double *resposta, Tempo_t *tempo, int *numIteracoes);
void FatorLU(SnlVar_t *var, int *trocas, int n);
void EliminacaoLU(SnlVar_t *var, int *trocas, int n);
void normsubs(SistLinear_t *SL, double *X);

#endif
