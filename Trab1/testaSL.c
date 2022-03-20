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

int main() {
    void *f;

    char *names[] = { "x1" };

    double X_new[] = { 0 }; 
    double X_old[1];
    double delta[1];
    int count = 1;

    f = evaluator_create("7*x1-log(x1)");
    assert(f);


    // F'(Xi) = J(Xi) ou  parciais numericas
    void *d = evaluator_derivative(f, "x1");

    // HESSIANA
    void *h = evaluator_derivative(d, "x1");

    X_old[0] = 1.0; // X0

    double Dxi1;
    int i;
    for( i = 0; i < 20; i++)
    {
        Dxi1 = evaluator_evaluate(d, 1, names, X_old);
        if( Dxi1 < 0.1)
            break;            

        // H(Xi)*(delta - Xi) = D(Xi)
        SistLinear_t *SL = alocaSistLinear(1);
        SL->n = 1;

        // calcula H(Xi)
        SL->A[0][0] = evaluator_evaluate(h, 1, names, X_old);
        // calcula D(Xi)
        SL->b[0] = Dxi1;
        
        // calcula w
        eliminacaoGauss(SL, delta);
        X_new[0] = X_old[0] + delta[0];

        printf("iter %d: Xnew %f  Xold %f   delt %f \n", i, X_new[0], X_old[0], delta[0]);

        if( delta[0] < 0.0000001)
            break;

        X_old[0] = X_new[0];
    }

    printf("delta:  %f\n", X_new[0]);


    // printf("f1(0.1) = %g\n", evaluator_evaluate( f, 1, names, X_new ) );
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