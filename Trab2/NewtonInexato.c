/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat 
********************************************************/

#include "NewtonInexato.h"

/*Codigo para efetuar o metodo de gauss seidel*/
void gauss_seidel(SistLinear_t *SL, double *X)
{
    int i,j,k,ite=0;
    double soma,xk,norma,diff,erro;
    erro=1e-6;
    norma=1.0+erro;

    /*Zera o vetor de Delta toda vez que a função começa*/
    for ( i = 0; i< SL->n; i++)
        X[i]=0;    
    
    /*Recalcula todos os Xs do sistema utilizando o metodo de gauss saidel 
    até a iteração maxima ou a norma ser menor que o erro*/
    for( k=0;norma>erro;++k)
    {
        norma=0.0;
        for (i = 0; i <SL->n; ++i)
        {
            for(soma=0,j=0;j<i;++j)
                soma+=SL->A[i][j]*X[j];
            for(j=i+1;j<SL->n;++j)
                soma+=SL->A[i][j]*X[j];
            xk=(SL->b[i]-soma)/SL->A[i][i];
            diff=fabs(xk-X[i]);
            if(diff>norma)
                norma=diff;
            X[i]=xk;
        }
        ite++;
        if(ite>50)
            break;
    }
}

void NewtonInexato(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter)
{    
    double tauxP, tauxder, tauxSL;
    int itr = 0;
    SnlVar_t *ni = alocaSnlVar(snl->chute, snl->n);
    // --------LOOP PRINCIPAL-------- //
    for(int i = 0; i < snl->iteracao; i++)
    {
        tauxP = timestamp();
        
        tauxder = timestamp();
        substituteX(snl, ni);                   // calcula H[X] e J[X]
        tauxder = timestamp() - tauxder;
        snl2sl(snl, ni);                        // copia dados de snl em sl
        
        resposta[i] = evaluator_evaluate(snl->f, snl->n, snl->names, ni->x0);
        
        tauxSL = timestamp();
        gauss_seidel(ni->sl,ni->delta);         // calcula H[X]*delta = - J[X]  // A*x = -b
        tauxSL = timestamp() - tauxSL;
        
        calcDelta(ni, snl->n);                  // X[i+1] = X[i] + delta[i]
        
        tauxP = timestamp() - tauxP;
        
        t->totalMetodo += tauxP;
        t->derivadas += tauxder;
        t->totalSL += tauxSL;
        itr++;

        if(Parada(snl, ni->delta))
            break;
    }
    *nIter = itr;
    liberaSnlVar(ni, snl->n);
}