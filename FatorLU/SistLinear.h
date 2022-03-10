#ifndef __SIST_LINEAR__
#define __SIST_LINEAR__

#include "utils.h"

/***************************
 ESTRUTURA SISTEMA LINEAR
****************************/
/*
  Método 3: vetor de ponteiros de linhas contíguas
  http://wiki.inf.ufpr.br/maziero/doku.php?id=prog2:alocacao_dinamica_de_matrizes&s[]=aloca%C3%A7%C3%A3o&s[]=de&s[]=matrizes
*/
typedef struct {
  unsigned int n;   // dimensão do SL
  double *M;        // vetor nxn de posições da matriz
  double **A;       // matriz dos coeficientes do SL (vetor de ponteiros para posições de M)
  double *b;        // termos independentes do SL
} SistLinear_t;

SistLinear_t *alocaSistLinear(unsigned int n);

void liberaSistLinear(SistLinear_t *SL);

// devolve um outro Sistema Linear que é uma cópia do SL
SistLinear_t *dupSL(SistLinear_t *SL);

// ordem de leitura: n A b
SistLinear_t *lerSistLinear();

// devolve um vetor com o resíduo do SL com incógnitas X
double *residuo(SistLinear_t *SL, double *X);

// faz o cálculo do resíduo usando variáveis do tipo long double
long double *long_residuo(SistLinear_t *SL, long double *X);

#endif
