/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
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

// #define MAXIT_REFINAMENTO 5

// #define DEBUG_FLAG


int main() {
    // testaSnL();
    SnlVar_t *np, *nm, *ni;
    SistNl_t *snl, *Psnl, *Msnl, *Isnl;
    double ptoPadrao, ptoModif, ptoInexato;

    double TtotalEG, TtotalLU, TtotalGS, 
    TderivadasEG, TderivadasLU, TderivadasGS, 
    TslEG, TslLU, TslGS;

    TtotalEG = TtotalLU = TtotalGS = 
    TderivadasEG = TderivadasLU = TderivadasGS = 
    TslEG = TslLU = TslGS = 0;

    while(snl = lerSistNL())
    {      
        genSistNaoLinear(snl);  // calcula Jacobiana e Hessiana
        printf("%i\n", snl->n);
        printf("%s\n", snl->funcao);

        // snl precisa de copia, muda o He

        // x0, x1, delta e SL para NEWTON PADRAO
        SnlVar_t *np = genSnlVar(snl);
        SistNl_t *Psnl = CopySnL(snl);

        // x0, x1, delta e SL para NEWTON MODIFICADO
        SnlVar_t *nm = genSnlVar(snl);
        SistNl_t *Msnl = CopySnL(snl);

        // x0, x1, delta e SL para NEWTON INEXATO
        SnlVar_t *ni = genSnlVar(snl);
        SistNl_t *Isnl = CopySnL(snl);


        printf("#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");

        // --------LOOP PRINCIPAL-------- //
        for(int i = 0; i < snl->iteracao; i++)
        {
            printf("%-12d \t| ", i); // imprime iteração

            // ELIMINACAO GAUSS / NEWTON PADRAO //
            if(fabs(minDelta(np->delta, snl->n)) >= Psnl->eps || i == 0){
                TtotalEG = timestamp();
                ptoPadrao = NewtonPadrao(Psnl, np);
                TtotalEG = timestamp() - TtotalEG;
            }


            /*
            // FATORACAO LU / NEWTON MODIFICADO //
            if(!(fabs(minDelta(nm->delta)) < Msnl->eps) || i == 0){
                TtotalLU = timestamp();
                ptoModif = NewtonModificado(Msnl, nm);
                TtotalLU = timestamp() - TtotalLU;
            }

            printCol(ptoModif, Msnl, nm);


            // GAUSS SEIDEL / NEWTON INEXATO //
            if(!(fabs(minDelta(ni->delta)) < Isnl->eps) || i == 0){
                TtotalGS = timestamp();
                ptoInexato = NewtonInexato(Isnl, ni);
                TtotalGS = timestamp() - TtotalGS;
            }

            printCol(ptoInexato, Isnl, ni);
            */


            printf("\n");
            // se max(normal(dos 3 deltas)) for < eps, break(TODO)
            if((fabs(minDelta(np->delta, snl->n)) < snl->eps) )
                break;
        }
        printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e  |\n", TtotalEG, TtotalLU, TtotalGS);
        printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e  |\n", TderivadasEG, TderivadasLU, TderivadasGS);
        printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e  |\n#\n\n", TslEG, TslLU, TslGS);

        liberaSnlVar(np);
        liberaSnlVar(nm);
        liberaSnlVar(ni);

        liberaSistNl(Psnl);
        liberaSistNl(Msnl);
        liberaSistNl(Isnl);
    }    
    return 0;
}


int omain()
{
    // testaSnL();

    SistNl_t *snl;
    snl = lerSistNL();
    genSistNaoLinear(snl);

    snlinfo(snl);
    // genNames(snl);

    // void *f = evaluator_create(snl->funcao);
    // assert(f);

    // genJacobiana(snl);
    // genHessiana(snl);

    printf("evaluator %f\n", evaluator_evaluate(snl->f, snl->n, snl->names, snl->chute));
    // printf("jacobiana %f\n", evaluator_evaluate(snl->Bf[0], 1, snl->names, snl->chute));
    // printf("hessiana  %f\n", evaluator_evaluate(snl->Hf[0][0], snl->n, snl->names, snl->chute));

    // funcao 1:
    // evaluator: 7.000
    // jacobiana: 6.000
    // hessiana:  1.000
}

