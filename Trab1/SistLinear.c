/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "SistLinear.h"

SistLinear_t *alocaSistLinear(unsigned int n) {

  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));

  if (SL) {
    SL->n = n;

    SL->A = (double **) malloc(sizeof(double *)*n);
    if (!(SL->A)) {
      free(SL);
      return NULL;
    }

    SL->M = (double *) malloc(sizeof(double)*n*n);
    if (!(SL->M)) {
      free(SL->A);
      free(SL);
      return NULL;
    }

    for (int i=0; i < n; ++i)
      SL->A[i] = SL->M + i*n;

    SL->b = (double *) malloc(sizeof(double)*n);
    if (!(SL->b)) {
      free(SL->M);
      free(SL->A);
      free(SL);
      return NULL;
    }


    /*ADICIONEI*/
    /////////////////////////////////
    SL->t = (int *) malloc(sizeof(int)*n); 
      if (!(SL->t)) {
        free(SL->M);
        free(SL->A);
        free(SL);
        return NULL;
      }            
      for(int i = 0; i < n; i++)       
        SL->t[i] = i;

    /* aloca vetor p/ trocas efetuadas em LU
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
    }

  return SL;
}

void liberaSistLinear(SistLinear_t *SL) {
  free(SL->b);
  free(SL->t);
  free(SL->M);
  free(SL->A);
  free(SL);
}

