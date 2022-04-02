/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./testaSNL < sistemas.dat
********************************************************/

#include "NewtonPadrao.h"

// ELIMINACAO GAUSS ou NEWTON PADRAO
void NewtonPadrao(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter)
{    
    // inicializa tudo que precisa
    // double *resposta = malloc(sizeof(double) * snl->iteracao);
    double tauxP, tauxder, tauxSL;

    int itr = 0;


    SnlVar_t *np = alocaSnlVar(snl->chute, snl->n);

    // --------LOOP PRINCIPAL-------- //
    for(int i = 0; i < snl->iteracao; i++)
    {
        
        if(!Parada(snl, np->delta)){
            tauxP = timestamp();
                    
            tauxder = timestamp();
            substituteX(snl, np);                   // calcula H[X] e J[X]
            tauxder = timestamp() - tauxder;
            snl2sl(snl, np);                        // copia dados de snl em sl

            resposta[i] = evaluator_evaluate(snl->f, snl->n, snl->names, np->x0);
            
            tauxSL = timestamp();
            eliminacaoGauss(np->sl, np->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
            tauxSL = timestamp() - tauxSL;

            // varinfo(*np, *snl);
            calcDelta(snl, np);                     // X[i+1] = X[i] + delta[i]

            // Devolve f(X), ponto critico estimado

            tauxP = timestamp() - tauxP;

            t->totalMetodo += tauxP;
            t->derivadas += tauxder;
            t->totalSL += tauxSL;
            itr++;
        }
        else
            break;
    }

    *nIter = itr;

    liberaSnlVar(np, snl->n);
    // return resposta;
}
