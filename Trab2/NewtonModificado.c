/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "NewtonModificado.h"


// FATORACAO LU ou NEWTON MODIFICADO
void NewtonModificado(SistNl_t *snl, double *resposta, Tempo_t *tempo, int *numIteracoes)
{
  SnlVar_t *nm = alocaSnlVar(snl->chute, snl->n);
  int *trocas = malloc(snl->n * sizeof(int));
  double tauxM, tauxder, tauxSL, aux;
  int itr = 0; 

  // --------LOOP PRINCIPAL-------- //
  for(int i = 0; i < snl->iteracao; i++)
  {
    tauxM = timestamp();

    if(i % snl->n == 0) 
    {
      tauxder = timestamp();
      calcHessiana(snl, nm);            // nm.He = snl.Hf[nm.X0]
      tauxder = timestamp() - tauxder;
      FatorLU(nm, trocas, snl->n);      // transforma nm. em LU
    }
    calcGradiente(snl, nm);             // calcula J[X]

    resposta[i] = evaluator_evaluate(snl->f, snl->n, snl->names, nm->x0); 

    tauxSL = timestamp();
    EliminacaoLU(nm, trocas, snl->n);     // resolve SL com LU
    tauxSL = timestamp() - tauxSL;
    calcDelta(nm, snl->n);                // X[i+1] = X[i] + delta[i]

    itr++;
    tauxM = timestamp() - tauxM;

    tempo->totalMetodo += tauxM;
    tempo->derivadas += tauxder;
    tempo->totalSL += tauxSL;

    if(Parada(snl, nm->delta) )
      break;    

  }
  *numIteracoes = itr;

  liberaSnlVar(nm, snl->n);
  free(trocas);
}

// calcula X com o sistema LUx = b
void EliminacaoLU(SnlVar_t *var, int *trocas, int n){
  
  // salva parametros em sistema linear temporario
  SistLinear_t *LU = alocaSistLinear(n);
  for(int i = 0; i < n; i++) for(int j = 0; j < n; j++)
    LU->A[i][j] = var->He[i][j];
  LU->n = n;

  // trocas no vetor (pivoteamento)
  double *aux = malloc(n * sizeof(double));
  for(int i = 0; i < n; i++)
    aux[i] = - var->Ge[ trocas[i] ];
  for(int i = 0; i < n; i++)
    LU->b[i] = aux[i];
  free(aux);

  // LZ = b
  double *Z = malloc(n * sizeof(double));
  normsubs(LU, Z); // LOWER (aplica trocas)

  // U*delta = Z 
  for (int j = 0; j < LU->n; ++j)
    LU->b[ j ] = Z[j];
  retrossubs(LU, var->delta);  // UPPER

  free(Z);
  liberaSistLinear(LU);
}


// TRIANGULACAO(triang) -> L = m  & U -> triangular normal
void FatorLU(SnlVar_t *var, int *trocas, int n){ 
  int max_i, aux;
  double m;

  // salva parametros em sistema linear temporario
  SistLinear_t *LU = alocaSistLinear(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      LU->A[i][j] = var->He[i][j];
    LU->b[i] = - var->Ge[i];
  }
  LU->n = n;
  
  for (int i = 0; i < LU->n; ++i) 
  {
    pivot(LU, i);

    for (int k = i+1; k < LU->n; ++k) 
    {
      m = LU->A[ k ][i] / LU->A[i][i];
      if (isnan(m))
        fprintf(stderr, "ERRO triang LU: %g ", LU->A[i][i]);
      // guarda m em L (ALTERA DEPOIS DA PRIMEIRA ITER)
      LU->A[ k ][i] = m;

      for (int j = i+1; j < LU->n; ++j)
        LU->A[ k ][j] -= LU->A[i][j] * m;
    }
  }

  for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++)
        var->He[i][j] = LU->A[i][j];
      trocas[i] = LU->t[i];
    }

  liberaSistLinear(LU);
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
