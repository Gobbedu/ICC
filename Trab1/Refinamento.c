#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"
#include "EliminacaoGauss.h"
#include "Refinamento.h"

// aplica MAXIT iterações do Refinamento em SL com solução X
void refinamento(SistLinear_t *SL, double *X, int MAXIT) {

  double *W = (double *) malloc(sizeof(double)*SL->n);

  long double *LX = (long double *) malloc(sizeof(long double)*SL->n); 
  for (int i = 0; i < SL->n; ++i)
    LX[i] = (long double) X[i];

  printf("  --> Refinamento:\n");
  for (int i = 0; i < MAXIT; ++i) {
    long double *res = long_residuo(SL, LX);
    printf("  ----> iteração %d. Resíduo: ", i);
    prnVetorLongDouble(res, SL->n);

    SistLinear_t *U = dupSL(SL);
    for (int j = 0; j < SL->n; ++j)
      U->b[j] = (double) res[j];

    eliminacaoGauss(U, W);

    for(int k = 0; k < SL->n; ++k)
      LX[k] += (long double) W[k];
    
    liberaSistLinear(U);
    free(res);
  }
}

void refinamentoLU(SistLinear_t *SL, double *X, int MAXIT) {

  SistLinear_t *LU;
  double *Z = (double *) malloc(sizeof(double)*SL->n);
  double *W = (double *) malloc(sizeof(double)*SL->n);

  long double *LX = (long double *) malloc(sizeof(long double)*SL->n); 
  for (int i = 0; i < SL->n; ++i)
    LX[i] = (long double) X[i];

  // chamar fatoracao LU aki,  PIVO HABILITADO
  LU = dupSL(SL);
  FatorLU(SL, LU);

  printf("  --> Refinamento LU:\n");
  for (int i = 0; i < MAXIT; ++i) 
  {
    long double *res = long_residuo(SL, LX);
    printf("  ----> iteração %d. Resíduo: ", i);
    prnVetorLongDouble(res, SL->n);

    // Aw = residuo
    // A = LU
    // LUw = r

    // L(z) = r // normsubs(L, z)
    // U(^x) = z , // retrossub(U, ^x)
    //onde ^x => w (Aw = r)

    //L(Z) = res
    for (int j = 0; j < LU->n; ++j)
      LU->b[ LU->t[j] ] = (double) res[j];
    normsubs(LU, Z);

    // UZ = w
    for (int j = 0; j < LU->n; ++j)
      LU->b[ LU->t[j] ] = (double) Z[j];
    retrossubs(LU, W);

    for(int k = 0; k < SL->n; ++k)
      LX[k] += (long double) W[k];
    
    free(res);
  }

  liberaSistLinear(LU);
}

// L e U presentes na msm matriz
void FatorLU(SistLinear_t *SL, SistLinear_t *LU)
{
  // TRIANGULACAO(triang) -> L = m  & U -> triangular normal
  int max_i, aux;
  double m;

  for (int i = 0; i < LU->n; ++i) 
  {
    pivot(LU, i);

    for (int k = i+1; k < LU->n; ++k) 
    {
      m = LU->A[k][i] / LU->A[i][i];
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
void normsubs(SistLinear_t *SL, double *X) {
  // para fatoracao LU, diagonal = 1
  for (int i = 0; i < SL->n; i++) {
    X[i] = SL->b[i];
    for (int j = 0; j < i; j++)
      X[i] -= SL->A[i][j] * X[j];
  }
}
