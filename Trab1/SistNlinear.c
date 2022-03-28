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

    // calcula f(X) antes de atualizar o valor
    // se valor existe, imprime na coluna
    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, np->x0), snl, np);
    
    snl2sl(snl, np->sl);                    // copia dados de snl em sl
    eliminacaoGauss(np->sl, np->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
    // varinfo(*np, *snl);
    calcDelta(np, snl->n);                  // X[i+1] = X[i] + delta[i]



    // printf("raiz funcao: ");
    // for(int i = 0; i < snl->n; i++)
    //     printf(" %g ", x1[i]);
    // printf("\n");

    // Devolve f(X), ponto critico estimado
    return evaluator_evaluate(snl->f, snl->n, snl->names, np->x1);
}

double NewtonModificado(SistNl_t *snl, SnlVar_t *nm)
{
    /* NEWTON MODIFICADO
 
    if(i % HESS_STEP == 0) 
        substituteX(snl, old_values); // calcula H[X] & J[X]

    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, np->x0), snl, np);

    snl2sl(snl, sl); // H[X]*delta = - J[X]  // A*x = -b
    // FATORACAO LU (TODO)
    eliminacaoGauss(sl, delta);

    calcDelta(new_values, old_values, delta, snl->n);
    */

    //PODE DELETAR AKI
    substituteX(snl, nm->x0);               // calcula H[X] e J[X]
    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, nm->x0), snl, nm);
    snl2sl(snl, nm->sl);                    // copia dados de snl em sl
    eliminacaoGauss(nm->sl, nm->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
    calcDelta(nm, snl->n);                  // X[i+1] = X[i] + delta[i]
    return evaluator_evaluate(snl->f, snl->n, snl->names, nm->x1);  
}

double NewtonInexato(SistNl_t *snl, SnlVar_t *nm)
{
    /* NEWTON INEXATO

    // GAUSS SEIDEL (TODO)

    */

    //PODE DELETAR AKI
    substituteX(snl, nm->x0);               // calcula H[X] e J[X]
    printCol(evaluator_evaluate(snl->f, snl->n, snl->names, nm->x0), snl, nm);
    snl2sl(snl, nm->sl);                    // copia dados de snl em sl
    eliminacaoGauss(nm->sl, nm->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
    calcDelta(nm, snl->n);                  // X[i+1] = X[i] + delta[i]
    return evaluator_evaluate(snl->f, snl->n, snl->names, nm->x1);  

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

void calcDelta(SnlVar_t *var, int n){
    for(int i = 0; i < n; i++){
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

void liberaSistNl(SistNl_t *snl) {
    for(int i = 0; i < snl->n; i++)
    {
        free(snl->names[i]);
        free(snl->He[i]);
    }

    free(snl->names);
    free(snl->He);
    free(snl->chute);
    free(snl->Be);
    
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

SnlVar_t *genSnlVar(SistNl_t *snl)
{
    SnlVar_t *var = malloc(sizeof(SnlVar_t));

    var->sl = alocaSistLinear(snl->n);
    var->x0 = genValues(snl->n, 0);
    var->x1 = genValues(snl->n, 0);
    var->delta = genValues(snl->n, 1);

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
    free(var);
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

double minDelta(double *delta, int n){
    double min = +INFINITY;

    for(int i = 0; i < n ; i++)
        min = (delta[i] < min) ? delta[i] : min;
        
    return min;
}

void snl2sl(SistNl_t *snl, SistLinear_t *sl)
{
    for(int i = 0; i < snl->n; ++i) {
        for(int j = 0; j < snl->n; ++j)
        {
            sl->A[i][j] = snl->He[i][j];
        }
        sl->b[i] = - snl->Be[i];
        /* dup->t[i] = SL->t[i]; // trocas ????? */
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
