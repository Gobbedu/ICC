#include "NewtonModificado.h"

// FATORACAO LU ou NEWTON MODIFICADO
void NewtonModificado(SistNl_t *snl, double *resposta, Tempo_t *tempo, int *numIteracoes)
{
  SnlVar_t *nm = alocaSnlVar(snl->chute, snl->n);
  double tauxM, tauxder, tauxSL, aux;
  int itr = 0; 

  // --------LOOP PRINCIPAL-------- //
  for(int i = 0; i < snl->iteracao; i++)
  {
    tauxM = timestamp();

    if(i % snl->n == 0) 
    {
      calcHessiana(snl, nm);    // calcula H[X]
      // substituteX(snl, nm);     // calcula H[X] e J[X]
      snl2sl(snl, nm);          // H[X]*delta = - J[X]  // A*x = -b
      FatorLU(nm->sl);          // transforma sl em LU
    }
    calcJacobiana(snl, nm);     // calcula J[X]

    resposta[i] = evaluator_evaluate(snl->f, snl->n, snl->names, nm->x0); 

    snl2sl(snl, nm);                    // H[X]*delta = - J[X]  // A*x = -b
    // varinfo(*nm, *snl);
    EliminacaoLU(nm->sl, nm->delta);    // resolve SL com LU
    // varinfo(*nm, *snl);

    calcDelta(snl, nm);                 // X[i+1] = X[i] + delta[i]



    itr++;
    tauxM = timestamp() - tauxM;
    tempo->totalMetodo += tauxM;

    if(Parada(snl, nm->delta) )
      break;    

  // if(i > 1 && fabs((resposta[i]-resposta[i-1])/resposta[i]) > 1 )  // se erro relativo > 1 PARE
  //   break;
  }
  *numIteracoes = itr;
  liberaSnlVar(nm, snl->n);
}


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
      m = LU->A[ k ][i] / LU->A[i][i];
      // if (isnan(m))
      //   printf("ERRO: %g ", LU->A[i][i]);
      // guarda m em L (ALTERA DEPOIS DA PRIMEIRA ITER)
      LU->A[ k ][i] = m;

      for (int j = i+1; j < LU->n; ++j)
        LU->A[ k ][j] -= LU->A[i][j] * m;
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
    normsubs(LU, Z); // LOWER (aplica trocas)

    // Ux = Z (as trocas ja foram feitas em Z)
    for (int j = 0; j < LU->n; ++j){
      LU->b[ j ] = (double) Z[j];
    }
    retrossubs(LU, X);  // UPPER

  free(Z);
}