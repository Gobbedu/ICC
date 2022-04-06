/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#ifndef __UTILS__
#define __UTILS__

#include <stdlib.h>
#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matheval.h>
#include <assert.h>

typedef struct {
    double totalSL; // tempo total do sistema linear & derivadas (matheval)
    double derivadas;
    double totalMetodo;
} Tempo_t;

void initTempo(Tempo_t *tempo);

double timestamp(void);

void prnVetorFloat(float *x, int n);

void prnVetorDouble(double *x, int n);

void prnVetorLongDouble(long double *x, int n);

#endif
