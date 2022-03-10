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
  double *W = (double *) malloc(sizeof(double)*SL->n);

  long double *LX, *Z = (long double *) malloc(sizeof(long double)*SL->n); 
  for (int i = 0; i < SL->n; ++i)
    LX[i] = (long double) X[i];

  // chamar fatoracao LU aki,  PIVO DESABILITADO
  FatorLU(SL, LU);

  printf("  --> Refinamento:\n");
  for (int i = 0; i < MAXIT; ++i) {
    long double *res = long_residuo(SL, LX);
    printf("  ----> iteração %d. Resíduo: ", i);
    prnVetorLongDouble(res, SL->n);

    // eliminacaoGauss(U, W);
    // L(z) = r // normsubs(L, z)
    // U(^x) = z , // retrossub(U, ^x)
    //onde ^x => w (Aw = r)

    //L(Z) 
    for (int j = 0; j < LU->n; ++j)
      LU->b[j] = (double) res[j];
    normsubs(LU, res);

    for (int j = 0; j < LU->n; ++j)
      Z[j] = (double) LU->b[j];

    // passa Z para LU->b
    retrossubs(LU, Z)
    // LU->b === W

    for(int k = 0; k < SL->n; ++k)
      LX[k] += (long double) LU->b[k];
    
    liberaSistLinear(LU);
    free(res);
  }
}

// L e U presentes na msm matriz
void FatorLU(SistLinear_t *SL, SistLinear_t *LU)
{
  LU = dupSL(SL);

  // TRIANGULACAO -> L = m  & U -> triangular normal
    for (int i = 0; i < LU->n; ++i) {
    // pivot(LU, i);
    for (int k = i+1; k < LU->n; ++k) {
      double m = LU->A[k][i] / LU->A[i][i];
      if (isnan(m))
        printf("ERRO: %g\n", LU->A[i][i]);
      // LU->A[k][i] = 0.0; // guarda m em L (Ld = 1)
      LU->A[k][i] = m;

      for (int j = i+1; j < LU->n; ++j)
        LU->A[k][j] -= LU->A[i][j] * m;
      LU->b[k] -= LU->b[i] * m;
    }
  }


}