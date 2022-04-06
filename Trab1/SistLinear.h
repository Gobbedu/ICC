/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

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
  int *t;           // trocas efetuadas em LU
} SistLinear_t;

SistLinear_t *alocaSistLinear(unsigned int n);

void liberaSistLinear(SistLinear_t *SL);

#endif
