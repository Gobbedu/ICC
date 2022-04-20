/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
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

#include "Rosenbrock.h"

// #define LIKWID_PERFMONI
#ifdef LIKWID_PERFMONI
#include <likwid.h>
#endif

typedef struct {
    double totalSL; // tempo total do sistema linear & derivadas (matheval)
    double derivadas;
    double Gradiente;
    double Hessiana;
    double totalMetodo;
} Tempo_t;

void initTempo(Tempo_t *tempo);

double timestamp(void);

void prnVetorFloat(float *x, int n);

void prnVetorDouble(double *x, int n);

void prnVetorLongDouble(long double *x, int n);

#endif
