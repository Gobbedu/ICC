#ifndef __UTILS__
#define __UTILS__

#include <stdlib.h>
#include <sys/time.h>

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


void copyDoubleVetor(double *entrada, double *destino, int n);

void copyDoubleMatrix(double **entrada, double **saida, int n);

void copyVoidMatrix(void ***entrada, void ***saida, int n);


#endif
