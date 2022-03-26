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

/* FUNCIONA */
void testaSnL(void){
    void *f;

    char *names[] = { "x1" };

    double X_new[] = { 0 }; 
    double X_old[1];
    double delta[1];
    double A;
    double b;

    int count = 1;

    f = evaluator_create("7*x1-log(x1)"); assert(f);

    // F'(Xi) = J(Xi) ou  parciais numericas
    void *d = evaluator_derivative(f, "x1");

    // HESSIANA
    void *h = evaluator_derivative(d, "x1");

    SistLinear_t *SL = alocaSistLinear(1);
    SL->n = 1;
    X_old[0] = 0.1;

    SL->A[0][0] = 1.0;

    int i;
    for( i = 0; i < 20; i++)
    {
        if( SL->A[0][0] < 0.1 )
            break;            

		// x_new = x_old - (Px/Dx); 
        // deriva tudo mais uma vez
        // H(Xi)*(delta - Xi) =  - D(Xi)
        
        SL->A[0][0] = evaluator_evaluate(h, 1, names, X_old);
        SL->b[0] = - (evaluator_evaluate(d, 1, names, X_old));
        // SL->b[0] = 0.01;
        
        // calcula w
        // SL->A * w = SL->b
        eliminacaoGauss(SL, delta);

        printf("A: %f   w: %f  b: %f\n", SL->A[0][0], delta[0], SL->b[0]);
        X_new[0] = X_old[0] + delta[0];

        printf("iter %d: Xnew %f  Xold %f   delt %f \n", i, X_new[0], X_old[0], delta[0]);

        X_old[0] = X_new[0];

        if( delta[0] < 0.0000001)
            break;
    }
    liberaSistLinear(SL);
    printf("raiz:  %f\n", X_new[0]);
}

void genNames(SistNl_t *snl){
    snl->names = malloc(snl->n * sizeof(char *));
    
    for(int i = 0; i < snl->n; i++){
        snl->names[i] = malloc(2 * sizeof(char));

        snl->names[i][0] = 'x';
        snl->names[i][1] = i+1+'0';
    }
}
// freeNames

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
