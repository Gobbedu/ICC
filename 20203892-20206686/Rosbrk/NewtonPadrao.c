/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "NewtonPadrao.h"

// // ELIMINACAO GAUSS ou NEWTON PADRAO
void NewtonPadrao(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter)
{    
	// inicializa tudo que precisa
	// double *resposta = malloc(sizeof(double) * iteracao);
	double tTotal, tauxder, tGrad, tHess, tauxSL;

	int n = snl->n;
	int itr = 0;

	char *mMetodo, *mGrad, *mHess, *mSL;
	mMetodo = malloc(sizeof(char)*50);
	mGrad   = malloc(sizeof(char)*50);
	mHess   = malloc(sizeof(char)*50);
	mSL     = malloc(sizeof(char)*50);

	sprintf(mMetodo, "PadraoMETODO_%u", n);
	sprintf(mGrad, "PadraoGRAD_%u", n);
	sprintf(mHess, "PadraoHESS_%u", n);
	sprintf(mSL, "PadraoSISTLIN_%u", n);

	SnlVar_t *np = alocaSnlVar(snl->chute, n);

	#ifdef LIKWID_PERFMONI
	LIKWID_MARKER_START(mMetodo);
	#endif

	// --------LOOP PRINCIPAL-------- //
	for(int i = 0; i < snl->iteracao; i++)
	{
		tTotal = timestamp();

		// tauxder = timestamp();
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mGrad);
		#endif
		tGrad = timestamp();
		calcGradiente(snl, np);					// calcula J[X]
		tGrad = timestamp() - tGrad;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mGrad);			
		#endif

		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mHess);
		#endif
		tHess = timestamp();
		calcHessiana(snl, np);					// calcula H[X]
		tHess = timestamp() - tHess;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mHess);			
		#endif
		// tauxder = timestamp() - tauxder;

		snl2sl(snl, np);                        // copia dados de snl em sl

		resposta[i] = rosenbrock(np->x0, n);
		
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mSL);
		#endif
		tauxSL = timestamp();
		eliminacaoGauss(np->sl, np->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
		tauxSL = timestamp() - tauxSL;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mSL);			
		#endif

		// varinfo(*np, *snl);
		calcDelta(np, n);                     // X[i+1] = X[i] + delta[i]

		// Devolve f(X), ponto critico estimado

		tTotal = timestamp() - tTotal;

		t->totalMetodo += tTotal;
		t->Gradiente += tGrad;
		t->Hessiana += tHess;
		t->totalSL += tauxSL;
		// t->derivadas += tauxder;
		itr++;
		
			if(Parada(snl, np->delta))
			break;
	}
	#ifdef LIKWID_PERFMONI
	LIKWID_MARKER_STOP(mMetodo);
	#endif

	free(mMetodo);
	free(mGrad);
	free(mHess);
	free(mSL);

	*nIter = itr;

	liberaSnlVar(np, n);
	// return resposta;
}


void pivot(SistLinear_t *SL, int i) {
	double max = fabs(SL->A[i][i]);
	int max_i = i;

	for (int j = i+1; j < SL->n; ++j) 
	{
		double v = fabs(SL->A[j][i]);
		if (v > max) {
			max = v;
			max_i = j;
		}
	}

	if (max_i != i) 
	{
		double *tmp = SL->A[i];
		SL->A[i] = SL->A[max_i];
		SL->A[max_i] = tmp;

		double aux = SL->b[i];
		SL->b[i] = SL->b[max_i];
		SL->b[max_i] = aux;

		int iaux = SL->t[i];
		SL->t[i] = SL->t[max_i];
		SL->t[max_i] = iaux;
	}
} 

// calcula retrossub com SL-b e salva em X
void retrossubs(SistLinear_t *SL, double *X) 
{
	for (int i = SL->n-1; i >=0; --i) 
	{
		X[i] = SL->b[i];
		for (int j = i+1; j < SL->n; j++)
			X[i] -= SL->A[i][j] * X[j];
		X[i] /= SL->A[i][i];
	}
}

void triang(SistLinear_t *SL) 
{
	for (int i = 0; i < SL->n; ++i) 
	{
		pivot(SL, i);
		for (int k = i+1; k < SL->n; ++k) 
		{
			double m = SL->A[k][i] / SL->A[i][i];
			if (isnan(m))
			  fprintf(stderr, "ERRO: %g ", SL->A[i][i]);
			SL->A[k][i] = 0.0;
			for (int j = i+1; j < SL->n; ++j)
				SL->A[k][j] -= SL->A[i][j] * m;
			SL->b[k] -= SL->b[i] * m;
		}
	}
}

void eliminacaoGauss(SistLinear_t *SL, double *X) 
{
	triang(SL);
	retrossubs(SL, X);
}

