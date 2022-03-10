#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"
#include "EliminacaoGauss.h"
#include "Refinamento.h"

// aplica MAXIT iterações do Refinamento em SL com solução X
void default_refinamento(SistLinear_t *SL, double *X, int MAXIT) {

  double *W = (double *) malloc(sizeof(double)*SL->n);

  // LX == ultima solucao X do SL
  long double *LX = (long double *) malloc(sizeof(long double)*SL->n); 
  for (int i = 0; i < SL->n; ++i)
    LX[i] = (long double) X[i];

  printf("  --> Refinamento:\n");
  for (int i = 0; i < MAXIT; ++i) {
    // residuo == res
    long double *res = long_residuo(SL, LX);
    printf("  ----> iteração %d. Resíduo: ", i);
    prnVetorLongDouble(res, SL->n);

    // sist U contem residuo em b, == Aw = r
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

void refinamento(SistLinear_t *SL, double *X, int MAXIT) {

  double *W = (double *) malloc(sizeof(double)*SL->n);

  // LX == ultima solucao X do SL
  long double *LX = (long double *) malloc(sizeof(long double)*SL->n); 
  for (int i = 0; i < SL->n; ++i)
    LX[i] = (long double) X[i];

  // A = LU
  // SL->A = L * U
  // fatorLU(SL, L, U, operL)

  printf("  --> Refinamento:\n");
  for (int i = 0; i < MAXIT; ++i) {

    /* residuo res = b - Ax
     * res = SL->b - SL->A * LX
     * RESIDUO:                 */
    long double *res = long_residuo(SL, LX);
    printf("  ----> iteração %d. Resíduo: ", i);
    prnVetorLongDouble(res, SL->n);


    // calcula w usando LU com == Aw = r
    SistLinear_t *E = dupSL(SL);
    for (int j = 0; j < SL->n; ++j)
      E->b[j] = (double) res[j];
    // A->b = b
    // E->b = res

    // eliminacaoGauss(E, W); // -> triang(SL); retrossub(E, W)
    
    // A = LU
    // LUw = r
    // Uw = z
    // Lz = r
    // fatorLU(inp A, )

    // LX == solucao x com ajuste W
    for(int k = 0; k < SL->n; ++k)
      LX[k] += (long double) W[k];
    
    liberaSistLinear(E);
    free(res);
  }
}

// void FatoracaoLU(SistLinear_t *SL, SistLinear_t *L, SistLinear_t *U, SistLinear_t *operL)

