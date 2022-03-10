/********************************************************
*    Resolução de Sistemas Lineares
*    profs. Armando Delgado e Guilherme Derenievicz
*    CI1164 - DInf/UFPR
*
*    make ou make LONGDOUBLE=1
*    ./testaSL < sistemas.dat
********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistLinear.h"
#include "EliminacaoGauss.h"
#include "Refinamento.h"

#define MAXIT_REFINAMENTO 5

int main() {
	SistLinear_t *SL;
  int i = 0;

  while (SL = lerSistLinear()) {
    double *X = (double *) malloc(sizeof(double)*SL->n);

    SistLinear_t *CP = dupSL(SL);

    printf("\n***** Sistema %d --> n = %d\n", i+1, SL->n);

    double tTotal = timestamp();
    eliminacaoGauss(SL, X);
    tTotal = timestamp() - tTotal;
    printf("  --> X: ");
    prnVetorDouble(X, SL->n);
    printf("  --> Tempo: %lf ms\n", tTotal);

    double *res = residuo(SL, X);
    printf("  --> Residuo: ");
    prnVetorDouble(res, CP->n);

    tTotal = timestamp();
    refinamento(CP, X, MAXIT_REFINAMENTO);
    tTotal = timestamp() - tTotal;
    printf("  --> Tempo: %lf ms\n", tTotal);

    free(X);
    free(res);
    liberaSistLinear(SL);
    liberaSistLinear(CP);
    ++i;
  }
}

