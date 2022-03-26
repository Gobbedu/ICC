#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <matheval.h>
#include <assert.h>

#include "utils.h"
#include "SistLinear.h"
#include "EliminacaoGauss.h"
#include "Refinamento.h"
#include "SistNlinear.h"

double NewtonPadrao(SistNl_t *snl, SnlVar_t *np)
{
    substituteX(snl, np->x0);               // calcula H[X] e J[X]
    snl2sl(snl, np->sl);                    // copia dados de snl em sl
    eliminacaoGauss(np->sl, np->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b

    calcDelta(np->x1, np->x0, np->delta, snl->n);

    // printf("raiz funcao: ");
    // for(int i = 0; i < snl->n; i++)
    //     printf(" %g ", x1[i]);
    // printf("\n");

    // Devolve f(X), ponto critico estimado
    return evaluator_evaluate(snl->f, snl->n, snl->names, np->x1);
}

double NewtonModificado()
{
    /* NEWTON MODIFICADO
 
    if(i % HESS_STEP == 0) 
        substituteX(snl, old_values); // calcula H[X] & J[X]

    snl2sl(snl, sl); // H[X]*delta = - J[X]  // A*x = -b
    // FATORACAO LU (TODO)
    eliminacaoGauss(sl, delta);

    calcDelta(new_values, old_values, delta, snl->n);


    */
   return 0;
}

double NewtonInexato()
{
    /* NEWTON INEXATO

    snl2sl(snl, sl); // H[X]*delta = - J[X]  // A*x = -b

    // GAUSS SEIDEL (TODO)

    calcDelta(new_values, old_values, delta, snl->n);

  
    */
   return 0;
}


SistNl_t *CopySnL(SistNl_t *snl)
{
    SistNl_t *new = alocaSistNl(snl->n);

    new->eps = snl->eps;
    new->f = snl->f;
    new->iteracao = snl->iteracao;
    new->n = snl->n;

    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            new->Hf[i][j] = snl->Hf[i][j];
            new->He[i][j] = snl->He[i][j];
        }
        new->funcao[i] = snl->funcao[i];
        new->chute[i] = snl->chute[i];
        new->names[i] = snl->names[i];
        new->Bf[i] = snl->Bf[i];
        new->Be[i] = snl->Be[i];
    }
    
    return new;
}

SnlVar_t *genSnlVar(SistNl_t *snl)
{
    SnlVar_t *var = malloc(sizeof(SnlVar_t));

    var->sl = alocaSistLinear(snl->n);
    var->x0 = genValues(snl->n, 0);
    var->x1 = genValues(snl->n, 0);
    var->delta = genValues(snl->n, 0);

    // chute inicial
    for(int i = 0; i < snl->n; i++) 
        var->x0[i] = snl->chute[i]; 

    return var;
}

void liberaSnlVar(SnlVar_t *var)
{
    liberaSistLinear(var->sl);
    free(var->x0);
    free(var->x1);
    free(var->delta);
}

void snl2sl(SistNl_t *snl, SistLinear_t *sl)
{
    for(int i = 0; i < snl->n; ++i) {
        for(int j = 0; j < snl->n; ++j)
            sl->A[i][j] = snl->He[i][j];
        sl->b[i] = - snl->Be[i];
        /* dup->t[i] = SL->t[i]; // trocas ????? */
    }
}

SistNl_t *alocaSistNl(unsigned int n){
    SistNl_t *SnL = (SistNl_t *) malloc(sizeof(SistNl_t));

    if (SnL){
        SnL->n = n;

        // tratar !SnL -> free
        // SnL->names = (char **) malloc(sizeof(char *) * n);
        // for(int i = 0; i < n; i++)
        //     SnL->names[i] = malloc( 2*sizeof(char));

        // vetor de chutes
        SnL->chute = malloc(sizeof(double)*n);
        SnL->names = malloc(SnL->n * sizeof(char *));

        SnL->Hf = (void ***) malloc(sizeof(void **)*n);
        if (!(SnL->Hf)) {
            free(SnL);
            return NULL;
        }
        for (int i=0; i < n; ++i)
            SnL->Hf[i] = (void **) malloc(sizeof(void *)*n);


        SnL->He = (double **) malloc(sizeof(double*)*n);
        for(int i = 0; i < n; i++)
            SnL->He[i] = malloc(sizeof(double)*n);
        if (!(SnL->He)) {
            free(SnL->Hf);
            free(SnL);
            return NULL;
        }

        SnL->Bf = (void **) malloc(sizeof(void *)*n);
        if(!(SnL->Bf)){
            free(SnL->Hf);
            free(SnL->He);
            free(SnL);
        }

        SnL->Be = (double *) malloc(sizeof(double)*n);
        if (!(SnL->Be)) {
            free(SnL->Hf);
            free(SnL->He);
            free(SnL);
            return NULL;
        }
    }
  return SnL;
}

void liberaSistNl(SistNl_t *SL) {
  free(SL->names);
  free(SL->Be);
  free(SL->Bf);
  free(SL->He);
  free(SL->Hf);
  free(SL);
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

void genNames(SistNl_t *snl){    
    for(int i = 0; i < snl->n; i++){
        snl->names[i] = malloc(2 * sizeof(char));

        snl->names[i][0] = 'x';
        snl->names[i][1] = i+1+'0';
    }
}

double *genValues(int n, double init){
    double *values = malloc(n * sizeof(double));

    for(int i = 0; i < n; i++)
        values[i] = init;

    return values;
}
// freeValues

void genJacobiana(SistNl_t *snl){
    char name[2];
    name[0] = 'x';

    for(int i = 0; i < snl->n; i++){
        name[1] = i+1+'0';
        // printf("derivating in %s\n", name);
        snl->Bf[i] = evaluator_derivative(snl->f, name);
    }
}

void genHessiana(SistNl_t *snl){
    char name[2];
    name[0] = 'x';

    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            // deriva Jacob[i] n vezes em x1..xn
            name[1] = j+1+'0';
            snl->Hf[i][j] = evaluator_derivative(snl->Bf[i], name);
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

void substituteX(SistNl_t *snl, double *X){

    // H[X]
    for(int i = 0; i < snl->n ; i++)
    {
        for(int j = 0; j < snl->n; j++)
            // snl->He[i][j] = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, X);
            snl->He[i][j] = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, X);

        // J[X]
        snl->Be[i] = evaluator_evaluate(snl->Bf[i], snl->n, snl->names, X);
    }

}

void calcDelta(double *new_values, double *old_values, double *delta, int n){
    for(int i = 0; i < n; i++){
        new_values[i] = old_values[i] + delta[i];
        old_values[i] = new_values[i];
    }
}

double minDelta(double *delta){
    double min = INFINITY;

    for(int i = 0; delta[i] ; i++)
        min = (delta[i] < min) ? delta[i] : min;
        
    return min;
}

void snlinfo(SistNl_t *S){
    printf("%i\n", S->n);
    printf("%s\n#\n", S->funcao);

    /*
    printf("\n-------------SNL INFO-------------\n");
    printf("#  n: %i\n", S->n);
    printf("#  funcao: %s\n", S->funcao);
    printf("#  chute:");

    for(int i = 0; i < S->n; i++)
        printf(" %g ", S->chute[i]);

    printf("#  eps: %g\n", S->eps);
    printf("#  iteracao: %i\n", S->iteracao);

    for(int i = 0; i < S->n ; i++)
        printf("#  name %i: %s\n", i, S->names[i]);

    // JACOBIANA
    printf("\nJACOBIANA\n#");
    for(int i = 0; i < S->n; i++)
        printf("  %f   ", S->Be[i]);
    // for(int i = 0 ; i < snl->n; i++)
    //     printf("jacob %i : %f\n",i , evaluator_evaluate(snl->Bf[i], snl->n, names, values));

        
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

    printf("---------------------------------\n");
    */
}

void printCol(double pto, SistNl_t *snl, SnlVar_t *nt)
{
    if (isnan(pto) || isinf(pto))
        printf("%1.14e\t\t\t| ", pto);
    else if (fabs(minDelta(nt->delta)) >= snl->eps)
        printf("%1.14e\t| ", pto);
    else
        printf("\t\t\t| ");
}