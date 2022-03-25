/********************************************************
*    Resolução de Sistemas Lineares
*    profs. Armando Delgado e Guilherme Derenievicz
*    CI1164 - DInf/UFPR
*
*    make ou make LONGDOUBLE=1
*    ./testaSL < sistemas.dat
********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <matheval.h>
#include <assert.h>

#include "utils.h"
#include "SistLinear.h"
#include "EliminacaoGauss.h"
#include "Refinamento.h"

#define MAXIT_REFINAMENTO 5
#define HESS_STEP 1
#define STR_BUFFER 500

typedef struct {
    int n;
    double eps1;
    double eps2;
    int iteracao;
    char funcao[STR_BUFFER];
    
    void ***Hf;
    double *He;

    void **Bf;
    double *Be;
}SistNl_t;

SistNl_t *alocaSistNl(unsigned int n){
    SistNl_t *SnL = (SistNl_t *) malloc(sizeof(SistNl_t));

    if (SnL){
        SnL->n = n;

        SnL->Hf = (void ***) malloc(sizeof(void **)*n);
        if (!(SnL->Hf)) {
            free(SnL);
            return NULL;
        }
        for (int i=0; i < n; ++i)
            SnL->Hf[i] = (void **) malloc(sizeof(void *)*n);


        SnL->He = (double *) malloc(sizeof(double)*n);
        if (!(SnL->He)) {
        free(SnL->Hf);
        free(SnL);
        return NULL;
        }

        SnL->Bf = (void *) malloc(sizeof(void)*n);
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

        scanf("%lf", &(SnL->eps1));
        scanf("%lf", &(SnL->eps2));
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

char ** genNames(int n){
    char **names = malloc(n * sizeof(char *));
    for(int i = 0; i < n; i++){
        names[i] = malloc(2 * sizeof(char));

        names[i][0] = 'x';
        names[i][1] = i+1+'0';
    }

    // for(int i = 0; i < n; i++){
    //     printf("name %i : %s\n", i, names[i]);
    // }

    return names;
}
double *genValues(int n, double init){
    double *values = malloc(n * sizeof(double));

    for(int i = 0; i < n; i++)
        values[i] = init;

    return values;
}

void genJacobiana(void *f, SistNl_t *snl){
    void **Jacobiana = snl->Bf;

    char name[2];
    name[0] = 'x';

    for(int i = 0; i < snl->n; i++){
        name[1] = i+1+'0';
        printf("derivating in %s\n", name);
        Jacobiana[i] = evaluator_derivative(f, name);
    }

    snl->Bf = Jacobiana;
}

void genHessiana(SistNl_t *snl){
    void ***Hessiana = snl->Hf;
    char name[2];
    name[0] = 'x';

    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            // deriva Jacob[i] n vezes em x1..xn
            name[1] = j+1+'0';
            Hessiana[i][j] = evaluator_derivative(snl->Bf[i], name);
        }
    }

    snl->Hf = Hessiana;
}

int main() {
    // testaSnL();

    // printf("%i\n", S->n);
    // printf("%s\n", S->funcao);
    // printf("%lf\n", S->eps1);
    // printf("%lf\n", S->eps2);
    // printf("%i\n", S->iteracao);

    //create
    //derivate (J)
    //derivate (H)

    SistNl_t *snl;
    snl = lerSistNL();

    void *f = evaluator_create(snl->funcao);
    double *values = genValues(snl->n, 1);
    char **names = genNames(snl->n);

    genJacobiana(f, snl);
    // for(int i = 0 ; i < snl->n; i++)
    //     printf("jacob %i : %f\n",i , evaluator_evaluate(snl->Bf[i], snl->n, names, values));
    
    genHessiana(snl);
    // printf("HESSIANA:\n");
    // for(int i = 0 ; i < snl->n; i++){
    //     for(int j = 0; j < snl->n; j++)
    //         printf("%3f  ", evaluator_evaluate(snl->Hf[i][j], snl->n, names, values));
    //     printf("\n");}


    // H[X]
    // J[X]
    // H[X]*delta = - J[X]
    
    return 0;
}

/*

// cabeçalho
printf("Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");

// para cada iteração
for (...) {
    printf("%d \t\t| ", i); // imprime iteração

    if (...) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(fx) || isinf(fx))
            printf("%1.14e\t\t\t| ", fx);
        else
            printf("%1.14e\t| ", fx);
    }
    else
        printf("\t\t\t| ");

    // repete para as outras duas colunas...
}

// imprimir os tempos
printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e\n", TtotalEG, TtotalLU, TtotalGS);
printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e\n", TderivadasEG, TderivadasLU, TderivadasGS);
printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e\n", TslEG, TslLU, TslGS);

*/