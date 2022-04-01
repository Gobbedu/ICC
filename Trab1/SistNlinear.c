#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <matheval.h>
#include <assert.h>

#include "utils.h"
#include "SistLinear.h"
#include "EliminacaoGauss.h"
#include "FatoracaoLU.h"
#include "SistNlinear.h"

double NewtonPadrao(SistNl_t *snl, SnlVar_t *np)
{
    substituteX(snl, np->x0);               // calcula H[X] e J[X]

    // calcula f(X) antes de atualizar o valor
    // se valor existe, imprime na coluna
    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, np->x0), snl, np);
    
    snl2sl(snl, np->sl);                    // copia dados de snl em sl
    eliminacaoGauss(np->sl, np->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
    // varinfo(*np, *snl);
    calcDelta(snl, np);                  // X[i+1] = X[i] + delta[i]

    // Devolve f(X), ponto critico estimado
    return evaluator_evaluate(snl->f, snl->n, snl->names, np->x1);
}

double NewtonModificado(SistNl_t *snl, SnlVar_t *nm, int i)
{
    // NEWTON MODIFICADO
 
    if(i % HESS_STEP == 0) 
    {
        substituteX(snl, nm->x0);   // calcula H[X] e J[X]

        // ERRO EM substituteX()
        // for(int i = 0; i < snl->n; i++){
        //     for(int j = 0; j < snl->n; j++){
        //         if( isnan(snl->He[i][j]) ){printf("\nTA AKI PORRA\n"); break;}
        //     }
        // }

        snl2sl(snl, nm->sl);        // H[X]*delta = - J[X]  // A*x = -b
        FatorLU(nm->sl);            // transforma sl em LU
    }

    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, nm->x0), snl, nm);

    snl2sl(snl, nm->sl);                // H[X]*delta = - J[X]  // A*x = -b
    EliminacaoLU(nm->sl, nm->delta);    // resolve SL com LU

    calcDelta(snl, nm);              // X[i+1] = X[i] + delta[i]

    return evaluator_evaluate(snl->f, snl->n, snl->names, nm->x1);  
}

/*Gauss seidel do Vods, n funciona 100%*/
void gauss_seidel(SistLinear_t *SL, double *X)
{
    int i,j,q,d;
  unsigned int n = SL->n;
  double *r =  malloc((SL->n) * sizeof (double)) ;
  double temp,sum,erroMaximo,erroCalculado;
  double **A = SL->A;
  double *b = SL->b;
  double *x = malloc((SL->n) * sizeof(double));
  //X = malloc((SL->n) * sizeof(double));

  // A[n][n] = Matriz principal (e->f)
  // b[n] = vetor_independente (e->termos_independentes)

  for(i=0;i<n;i++){
      r[i] = 0;
  }
  
  q = 0;
  do{
      erroCalculado = 0;
      q++;
      for(i=0;i<n;i++){
          sum = 0;
          for(j=0;j<n;j++){
              if(i != j){
                  sum = sum + (A[i][j] * r[j]);
              }
            }
          temp = (-1.0 / A[i][i]) * sum + b[i] / A[i][i];
          erroMaximo = fabs(temp - r[i]);
          r[i] = temp;
          if(erroMaximo > erroCalculado)
              erroCalculado = erroMaximo;
        }
  }while(erroCalculado >= 1e-6  && q<=50);

  for(i=0;i<n;i++){ //copiar dados calculados para *x
    x[i] = r[i];
  }

  //X= x; //copiar dados para estrutura

  free(r);
  free(x);
  //free(X);


}
double NewtonInexato(SistNl_t *snl, SnlVar_t *ni)
{
    /* NEWTON INEXATO

    // GAUSS SEIDEL (TODO)

    */

    //PODE DELETAR AKI
    substituteX(snl, ni->x0);               // calcula H[X] e J[X]
    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, ni->x0), snl, ni);
    snl2sl(snl, ni->sl);                    // copia dados de snl em sl
    gauss_seidel(ni->sl,ni->delta);         // calcula H[X]*delta = - J[X]  // A*x = -b
    calcDelta(snl, ni);                  // X[i+1] = X[i] + delta[i]
    return evaluator_evaluate(snl->f, snl->n, snl->names, ni->x1);  

}


SistNl_t *lerSistNL(void)
{
    unsigned int n;

    SistNl_t *SnL = NULL;
  
    if (scanf("%d",&n) != EOF) 
    {  
        SnL = alocaSistNl(n);
        if (!SnL) return NULL;
        
        scanf("%s", SnL->funcao);

        for(int i = 0; i < n; i++)
            scanf("%lf", &(SnL->chute[i]));

        scanf("%lf", &(SnL->eps));
        scanf("%i", &(SnL->iteracao));

    }
  
  return SnL;
}

int Parada(SistNl_t *snl, SnlVar_t *nt)
{
    // normal -> || J(X) || = max{ |Ji(X)|, 1 ≤ i ≤ n}
    double maxF = -INFINITY;
    double Ji;
    for(int i = 0; i < snl->n; i++)
    {
        Ji = evaluator_evaluate(snl->Bf[i], snl->n, snl->names, nt->x0);
        maxF = (fabs(Ji) > maxF) ? fabs(Ji) : maxF;
    }

    // normal -> || X || = max{ |xi|, 1 ≤ i ≤ n}
    double maxD = -INFINITY;

    for(int i = 0; i < snl->n; i++){
        maxD = (fabs(nt->delta[i]) > maxD) ? fabs(nt->delta[i]) : maxD;
    }

    // se ||delta|| < eps PARE
    // se ||Ji(X)|| < eps PARE
    // return (maxF < snl->eps) || (maxD < snl->eps);
    return (maxD < snl->eps);

}

void calcDelta(SistNl_t *snl, SnlVar_t *var){
    if(fabs(minDelta(var->delta, snl->n)) >= snl->eps)
    for(int i = 0; i < snl->n; i++){
        var->x1[i] = var->x0[i] + var->delta[i];
        var->x0[i] = var->x1[i];
    }
}

void genNames(SistNl_t *snl){  
    for(int i = 0; i < snl->n; i++){
        sprintf(snl->names[i], "x%i\0", i+1);
    }
}

double *genValues(int n, double init){
    double *values = malloc(n * sizeof(double));

    for(int i = 0; i < n; i++)
        values[i] = init;

    return values;
}

void genJacobiana(SistNl_t *snl){
    for(int i = 0; i < snl->n; i++){
        snl->Bf[i] = evaluator_derivative(snl->f, snl->names[i]);
    }
}

void genHessiana(SistNl_t *snl){
    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            // deriva Jacob[i] n vezes em x1..xn
            snl->Hf[i][j] = evaluator_derivative(snl->Bf[i], snl->names[j]);
        }
    }
}

void genSistNaoLinear(SistNl_t *snl)
{
        // precisa ser chamado nesta ordem
        snl->f = evaluator_create(snl->funcao);
        assert(snl->f);

        genNames(snl);          // nomes das variaveis x1, x2 .. xn
        genJacobiana(snl);   // derivadas de f c/ respeito a x1, x2 .. xn
        genHessiana(snl);       // possiveis combinacoes de segunda derivada

}

SistNl_t *alocaSistNl(unsigned int n){
    SistNl_t *SnL = (SistNl_t *) malloc(sizeof(SistNl_t));

    if (SnL){
        SnL->n = n;

        // vetor de chutes
        SnL->chute = malloc(sizeof(double)*n);

        // snl->names[i] = malloc(4 * sizeof(char)); // maximo 999 variaveis
        SnL->names = malloc(SnL->n * sizeof(char *));
        for(int i = 0; i < SnL->n; i++)
            SnL->names[i] = malloc(5 * sizeof(char));

        SnL->Hf = (void ***) malloc(sizeof(void **)*n);
        if (!(SnL->Hf)) {
            free(SnL);
            return NULL;
        }
        for (int i=0; i < n; ++i)
            SnL->Hf[i] = (void **) malloc(sizeof(void *)*n);


        SnL->Bf = (void **) malloc(sizeof(void *)*n);
        if(!(SnL->Bf)){
            free(SnL->Hf);
            free(SnL);
        }

    }
  return SnL;
}

void liberaSistNl(SistNl_t *snl) {
    for(int i = 0; i < snl->n; i++)
    {
        free(snl->names[i]);
    }

    free(snl->names);
    free(snl->chute);
    free(snl);
}

void liberaMatheval(SistNl_t *snl)
{
    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            evaluator_destroy(snl->Hf[i][j]);
        }
        evaluator_destroy(snl->Bf[i]);
        free(snl->Hf[i]);
    }
    evaluator_destroy(snl->f);
    free(snl->Bf);
    free(snl->Hf);
}

SnlVar_t *alocaSnlVar(double *chute, int n)
{
    SnlVar_t *var = malloc(sizeof(SnlVar_t));

    var->sl = alocaSistLinear(n);
    var->x0 = genValues(n, 0);
    var->x1 = genValues(n, 0);
    var->delta = genValues(n, 1);

    // aloca He e libera var caso erro
    var->He = (double **) malloc(sizeof(double*)*n);
    if (!(var->He)) {
        liberaSistLinear(var->sl);
        free(var->x0);
        free(var->x1);
        free(var->delta);
        free(var);
        fprintf(stderr, "erro ao alocar Hessiana exata em variaveis do sist_n_linear\n");
        return NULL;
    }
    for(int i = 0; i < n; i++){
        var->He[i] = malloc(sizeof(double)*n);
    }
        
    var->Je = (double *) malloc(sizeof(double)*n);
    if (!(var->Je)) {
        liberaSistLinear(var->sl);
        free(var->x0);
        free(var->x1);
        free(var->delta);
        free(var->He);
        free(var);
        fprintf(stderr, "erro ao alocar Jacobiana exata em variaveis do sist_n_linear\n");
        return NULL;
    }

    // chute inicial
    for(int i = 0; i < n; i++) 
        var->x0[i] = chute[i]; 

    return var;
}

void liberaSnlVar(SnlVar_t *var)
{
    liberaSistLinear(var->sl);
    free(var->x0);
    free(var->x1);
    free(var->delta);

    for(int i = 0; i < var->n; i++)
    {
        free(var->He[i]);
    }
    free(var->He);
    free(var->Je);

    free(var);
}

void substituteX(SistNl_t *snl, double *X){
    double aux;
    // H[X]
    for(int i = 0; i < snl->n ; i++)
    {
        for(int j = 0; j < snl->n; j++){
            aux = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, X);
            if(isnan(aux)) printf("TA AKI O ERROOOOOOOOOO\n");
            snl->He[i][j] = aux;
        }

        // J[X]
        snl->Be[i] = evaluator_evaluate(snl->Bf[i], snl->n, snl->names, X);
    }

}

void snl2sl(SistNl_t *snl, SistLinear_t *sl)
{
    for(int i = 0; i < snl->n; ++i) {
        for(int j = 0; j < snl->n; ++j)
        {
            sl->A[i][j] = snl->He[i][j];
        }
        sl->b[i] = - snl->Be[i];
        // sl->t[i] ; // trocas salvas no pivo, dentro de nt->t[] (nao altera)
    }
}

void snlinfo(SistNl_t *S){
    printf("\n------------------------------------SNL INFO------------------------------------\n");
    printf("#  n: %i\n", S->n);
    printf("#  funcao: %s\n", S->funcao);

    printf("#  chute:");
    for(int i = 0; i < S->n; i++)
        printf(" %g ", S->chute[i]);
    printf("\n");

    printf("#  eps: %g\n", S->eps);
    printf("#  iteracao: %i\n", S->iteracao);

    printf("#  names: ");
    for(int i = 0; i < S->n ; i++)
        printf(" %s ", S->names[i]);
    printf("\n");

    printf("\nVETOR CHUTE's:\n#");
    for(int i = 0; i < S->n; i++)
        printf("  %f, ", S->chute[i]);
    printf("\nEVALUATOR\n#  f(CHUTE) = %g\n",evaluator_evaluate(S->f, S->n, S->names, S->chute));

    // JACOBIANA
    printf("\nJACOBIANA\n#");
    for(int i = 0; i < S->n; i++)
        printf(" %f ", S->Be[i]);
    printf("\n");
    // for(int i = 0 ; i < snl->n; i++)
    //     printf("jacob %i : %f\n",i , evaluator_evaluate(snl->Bf[i], snl->n, names, values));

        
    /*
    // HESSIANA
    printf("\nHESSIANA\n#");
    for(int i = 0; i < S->n; i++){
        for(int j = 0; j < S->n; j++){
            printf("  %f  ", S->He[i][j]);
        }
        printf("\n");
    }
    // for(int i = 0 ; i < snl->n; i++){
    //     for(int j = 0; j < snl->n; j++)
    //         printf("%3f  ", evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, old_values));
    //     printf("\n");}
    */

    printf("--------------------------------------END SNL INFO-----------------------------------------\n");
}

void varinfo(SnlVar_t nt, SistNl_t snl)
{
    printf("\n------------------------------------INFO VARIAVEIS------------------------------------\n");
    printf("x0: "); for(int i = 0; i < snl.n; i++) printf(" %f ", nt.x0[i]); printf("\n");
    printf("x1: "); for(int i = 0; i < snl.n; i++) printf(" %f ", nt.x1[i]); printf("\n");
    printf("delta: "); for(int i = 0; i < snl.n; i++) printf(" %f ", nt.delta[i]); printf("\n");
    printf("evaluator[x0] = %1.14e\n", evaluator_evaluate(snl.f, snl.n, snl.names, nt.x0));
    printf("-------------------------------------------------------------------------------\n");
}

void printCol(double pto, SistNl_t *snl, SnlVar_t *nt)
{
    if (isnan(pto) || isinf(pto))
        printf("%1.14e\t\t\t| ", pto);
    else if (fabs(minDelta(nt->delta, snl->n)) >= snl->eps)
        printf("%1.14e\t| ", pto);
    else
        printf("\t\t\t| ");
}
