#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef __SIST_LINEAR__
#include "SistLinear.h"
#endif

#ifndef __ELIM_GAUSS__
#include "EliminacaoGauss.h"
#endif

#include "FatoracaoLU.h"

// L e U presentes na msm matriz
void FatorLU(SistLinear_t *LU){ 
  // TRIANGULACAO(triang) -> L = m  & U -> triangular normal
  int max_i, aux;
  double m;

  for (int i = 0; i < LU->n; ++i) 
  {
    pivot(LU, i);

    for (int k = i+1; k < LU->n; ++k) 
    {
      if(LU->A[k][i] != 0 && LU->A[i][i] != 0)
        m = LU->A[k][i] / LU->A[i][i];
      else m = 0;
      if (isnan(m))
        printf("ERRO: %g ", LU->A[i][i]);
      // guarda m em L (ALTERA DEPOIS DA PRIMEIRA ITER)
      LU->A[ k ][i] = m;     

      for (int j = i+1; j < LU->n; ++j)
        LU->A[k][j] -= LU->A[i][j] * m;
      LU->b[k] -= LU->b[i] * m;
    }
  }

}

// calcula retrossub em L, com SL-b e salva em X
void normsubs(SistLinear_t *SL, double *X)  {
  // para fatoracao LU, diagonal = 1
  for (int i = 0; i < SL->n; i++) {
    X[i] = SL->b[i];
    for (int j = 0; j < i; j++)
      X[i] -= SL->A[i][j] * X[j];
  }
}

// calcula X com o sistema LUx = b
void EliminacaoLU(SistLinear_t *LU, double *X){
  double *Z = malloc(sizeof(double)* LU->n);

  // considerar trocas do pivoteamento
  //trocas nas posicoes de b 
    for (int i = 0; i < LU->n; ++i){
      LU->b[ LU->t[i] ] = LU->b[i];
    }
    normsubs(LU, Z);

    // Ux = Z (as trocas ja foram feitas em Z)
    for (int j = 0; j < LU->n; ++j){
      LU->b[ LU->t[j] ] = (double) Z[j];
    }
    retrossubs(LU, X);

  free(Z);
}