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
    double tTotal, tGrad, tHess, tauxder, tauxSL;
    int itr = 0;

    char *mMetodo, *mGrad, *mHess, *mSL;
	mMetodo = malloc(sizeof(char)*(20 + (log10(snl->n)+1)));
	mGrad   = malloc(sizeof(char)*(20 + (log10(snl->n)+1)));
	mHess   = malloc(sizeof(char)*(20 + (log10(snl->n)+1)));
	mSL     = malloc(sizeof(char)*(20 + (log10(snl->n)+1)));

	sprintf(mMetodo, "InexatoMETODO_%u", snl->n);
	sprintf(mGrad, "InexatoGRAD_%u", snl->n);
	sprintf(mHess, "InexatoHESS_%u", snl->n);
	sprintf(mSL, "InexatoSISTLIN_%u", snl->n);

    #ifdef LIKWID_PERFMONI
    LIKWID_MARKER_START(mMetodo);
    #endif
    SnlVar_t *ni = alocaSnlVar(snl->chute, snl->n);
    // --------LOOP PRINCIPAL-------- //
    for(int i = 0; i < snl->iteracao; i++)
    {
        tTotal = timestamp();
        
        // tauxder = timestamp();
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mHess);
		#endif
		tGrad = timestamp();
		calcGradiente(snl, ni);					// calcula J[X]
		tGrad = timestamp() - tGrad;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mHess);			
		#endif
        
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mGrad);
		#endif
		tHess = timestamp();
		calcHessiana(snl, ni);					// calcula H[X]
		tHess = timestamp() - tHess;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mGrad);			
		#endif
		// tauxder = timestamp() - tauxder;

        snl2sl(snl, ni);                        // copia dados de snl em sl
        
        resposta[i] = rosenbrock(ni->x0, snl->n);

		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mSL);
		#endif
        tauxSL = timestamp();
        gauss_seidel(ni->sl,ni->delta);         // calcula H[X]*delta = - J[X]  // A*x = -b
        tauxSL = timestamp() - tauxSL;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mSL);			
		#endif
        
        calcDelta(ni, snl->n);                  // X[i+1] = X[i] + delta[i]
        
        tTotal = timestamp() - tTotal;
        
        t->totalMetodo += tTotal;
		t->Gradiente += tGrad;
		t->Hessiana += tHess;
		t->totalSL += tauxSL;
		// t->derivadas += tauxder;
        itr++;

        if(Parada(snl, ni->delta))
            break;
    }

    #ifdef LIKWID_PERFMONI
    LIKWID_MARKER_STOP(mMetodo);
    #endif

    *nIter = itr;
    liberaSnlVar(ni, snl->n);
}