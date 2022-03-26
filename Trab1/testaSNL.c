/********************************************************
*    Resolução de Sistemas Lineares
*    profs. Armando Delgado e Guilherme Derenievicz
*    CI1164 - DInf/UFPR
*
*    make ou make LONGDOUBLE=1
*    ./testaSNL < sistemas.dat
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

#include "SistNlinear.h"

#define MAXIT_REFINAMENTO 5

// #define DEBUG_FLAG

void snl2sl(SistNl_t *snl, SistLinear_t *sl)
{
    for(int i = 0; i < snl->n; ++i) {
        for(int j = 0; j < snl->n; ++j)
            sl->A[i][j] = snl->He[i][j];
        sl->b[i] = - snl->Be[i];
        /* dup->t[i] = SL->t[i]; // trocas ????? */
    }
}

// NEWTON PADRAO
double NewtonPadrao(SistNl_t *snl)
{
    SistLinear_t *sl = alocaSistLinear(snl->n);

    double *x0 = genValues(snl->n, 0);
    double *x1 = genValues(snl->n, 0);
    double *delta = genValues(snl->n, 0);

    // chute inicial
    for(int i = 0; i < snl->n; i++) 
        x0[i] = snl->chute[i];

    for(int i = 0; i < snl->iteracao; i++)
    {
        substituteX(snl, x0);   // calcula H[X] e J[X]

        snl2sl(snl, sl);                // copia dados de snl em sl
        eliminacaoGauss(sl, delta);     // calcula H[X]*delta = - J[X]  // A*x = -b

        calcDelta(x1, x0, delta, snl->n);

        // calcular normal de delta (TODO)
        if(fabs(minDelta(delta)) < snl->eps)
            break;
    }
    liberaSistLinear(sl);

    // printf("raiz funcao: ");
    // for(int i = 0; i < snl->n; i++)
    //     printf(" %g ", x1[i]);
    // printf("\n");

    return evaluator_evaluate(snl->f, snl->n, snl->names, x1);
}
void NewtonModificado()
{
    /* NEWTON MODIFICADO
    for(int i = 0; i < snl->iteracao; i++)
    {
        if(i % HESS_STEP == 0) 
            substituteX(snl, old_values); // calcula H[X] & J[X]

        snl2sl(snl, sl); // H[X]*delta = - J[X]  // A*x = -b
        // FATORACAO LU (TODO)
        eliminacaoGauss(sl, delta);

        calcDelta(new_values, old_values, delta, snl->n);

        // calcular normal de delta (TODO)
        if(fabs(minDelta(delta)) < snl->eps)
            break;
    }*/

}
void NewtonInexato()
{
    /* NEWTON INEXATO
    for(int i = 0; i < snl->iteracao; i++)
    {

        snl2sl(snl, sl); // H[X]*delta = - J[X]  // A*x = -b
        // GAUSS SEIDEL (TODO)

        calcDelta(new_values, old_values, delta, snl->n);

        // calcular normal de delta (TODO)
        if(fabs(minDelta(delta)) < snl->eps)
            break;
    }*/

}

int main() {
    // testaSnL();
    SistNl_t *snl;
    while(snl = lerSistNL())
    {      
        genSistNaoLinear(snl);   // calcula Jacobiana e Hessiana
        snlinfo(snl);           // imprime dados do sistema

        printf("ponto critico:  %1.14e\n#\n", NewtonPadrao(snl));
    }    
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

int Testee()
{
    // testaSnL();

    SistNl_t *snl;
    snl = lerSistNL();
    double *valores = genValues(snl->n, 1);
    genSistNaoLinear(snl);
    // genNames(snl);

    // void *f = evaluator_create(snl->funcao);
    // assert(f);

    // genJacobiana(snl);
    // genHessiana(snl);

    printf("evaluator %f\n", evaluator_evaluate(snl->f, snl->n, snl->names, valores));
    printf("jacobiana %f\n", evaluator_evaluate(snl->Bf[0], 1, snl->names, valores));
    printf("hessiana  %f\n", evaluator_evaluate(snl->Hf[0][0], snl->n, snl->names, valores));
    // evaluator: 7.000
    // jacobiana: 6.000
    // hessiana:  1.000
}

