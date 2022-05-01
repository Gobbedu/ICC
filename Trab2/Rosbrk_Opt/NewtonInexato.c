/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat 
********************************************************/

#include "NewtonInexato.h"


void NewtonInexato(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter)
{    
    double tTotal, tGrad, tHess, tauxder, tauxSL;
    int n = snl->n;
    int itr = 0;


    char mMetodo[MARKER_SIZE], mGrad[MARKER_SIZE], 
	mHess[MARKER_SIZE], mSL[MARKER_SIZE];

	sprintf(mMetodo, "InexatoMETODO_%u", n);
	sprintf(mGrad, "InexatoGRAD_%u", n);
	sprintf(mHess, "InexatoHESS_%u", n);
	sprintf(mSL, "InexatoSISTLIN_%u", n);

    #ifdef LIKWID_PERFMONI
    LIKWID_MARKER_START(mMetodo);
    #endif
    SnlVar_t *ni = alocaSnlVar(snl->chute, n);
    // --------LOOP PRINCIPAL-------- //
    for(int i = 0; i < snl->iteracao; i++)
    {
        tTotal = timestamp();
        
        // tauxder = timestamp();
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mHess);
		#endif
		tGrad = timestamp();
        // calcula G[X]
        calcGradiente(ni->Ge, n, ni->x0);tGrad = timestamp() - tGrad;  		
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mHess);			
		#endif
        
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mGrad);
		#endif
		tHess = timestamp();
        // calcula H[X]
		calcHessiana(ni->He, n, ni->x0);tHess = timestamp() - tHess;   		
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mGrad);			
		#endif
		// tauxder = timestamp() - tauxder;

        resposta[i] = rosenbrock(ni->x0, n);

		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mSL);
		#endif
        tauxSL = timestamp();
        // calcula H[X]*delta = - J[X]  // A*x = -b
        gauss_seidel(ni, n); tauxSL = timestamp() - tauxSL;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mSL);			
		#endif
        
        calcDelta(ni, n);                  // X[i+1] = X[i] + delta[i]
        // X[i+1] = X[i] + delta[i]
		// for(int i = 0; i < n; i++)	ni->x1[i] = ni->x0[i] + ni->delta[i];
		// for(int i = 0; i < n; i++)	ni->x0[i] = ni->x1[i];

        
        tTotal = timestamp() - tTotal;
        
        t->totalMetodo += tTotal;
		t->Gradiente += tGrad;
		t->Hessiana += tHess;
		t->totalSL += tauxSL;
		// t->derivadas += tauxder;
        itr++;

        if(Parada(ni->delta, snl->eps, snl->n))
            break;
    }

    #ifdef LIKWID_PERFMONI
    LIKWID_MARKER_STOP(mMetodo);
    #endif

    *nIter = itr;
    liberaSnlVar(ni, n);
}


/*Codigo para efetuar o metodo de gauss seidel*/
void gauss_seidel(SnlVar_t *var, int n)
{
    int i,j,k,ite=0;
    double soma,xk,norma,diff,erro;
    erro=1e-6;
    norma=1.0+erro;

    /*Zera o vetor de Delta toda vez que a função começa*/
    for ( i = 0; i< n; i++)
        var->delta[i] = 0;    
    
    /*Recalcula todos os Xs do sistema utilizando o metodo de gauss saidel 
    até a iteração maxima ou a norma ser menor que o erro*/
    for( k=0;norma>erro;++k)
    {
        norma=0.0;
        for (i = 0; i < n; ++i)
        {
            for(soma=0,j=0;j<i;++j)
                soma+=var->He[i][j]*var->delta[j];

            for(j=i+1; j< n;++j)
                soma+=var->He[i][j]*var->delta[j];

            xk=(-var->Ge[i]-soma)/var->He[i][i];

            diff=fabs(xk-var->delta[i]);
            if(diff>norma)
                norma=diff;
            var->delta[i]=xk;
        }
        ite++;
        if(ite>50)
            break;
    }
}

