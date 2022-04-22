/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "NewtonPadrao.h"

// // ELIMINACAO GAUSS ou NEWTON PADRAO
// Devolve f(X), ponto critico estimado
void NewtonPadrao(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter)
{    
	double tTotal, tauxder, tGrad, tHess, tauxSL;
	int itr = 0;
	int n = snl->n;

	char mMetodo[MARKER_SIZE], mGrad[MARKER_SIZE], 
	mHess[MARKER_SIZE], mSL[MARKER_SIZE];

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
		// calcula G[X]
		calcGradiente(np->Ge, n, np->x0);tGrad = timestamp() - tGrad;					
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mGrad);			
		#endif

		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mHess);
		#endif
		tHess = timestamp();
		// calcula H[X]
		calcHessiana(np->He, n, np->x0);tHess = timestamp() - tHess;					
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mHess);			
		#endif
		// tauxder = timestamp() - tauxder;
		
		resposta[i] = rosenbrock(np->x0, n);
		
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_START(mSL);
		#endif
		tauxSL = timestamp();
		// calcula H[X]*delta = - J[X]  // A*x = -b
		eliminacaoGauss(np, n);     
		tauxSL = timestamp() - tauxSL;
		#ifdef LIKWID_PERFMONI
		LIKWID_MARKER_STOP(mSL);			
		#endif

		// X[i+1] = X[i] + delta[i]
		for(int i = 0; i < n; i++)	np->x1[i] = np->x0[i] + np->delta[i];
		for(int i = 0; i < n; i++)	np->x0[i] = np->x1[i];

		tTotal = timestamp() - tTotal;

		t->totalMetodo += tTotal;
		t->Gradiente += tGrad;
		t->Hessiana += tHess;
		t->totalSL += tauxSL;
		// t->derivadas += tauxder;
		itr++;
		
		if(Parada(np->delta, snl->eps, snl->n))
			break;
	}
	#ifdef LIKWID_PERFMONI
	LIKWID_MARKER_STOP(mMetodo);
	#endif

	*nIter = itr;

	liberaSnlVar(np, n);
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
	for (int i = SL->n-1; i >=0; --i) 	{
		X[i] = - SL->b[i];

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

void eliminacaoGauss(SnlVar_t *var, int n) 
{
	SistLinear_t SL;

	SL.A = var->He;
	SL.b = var->Ge;
	SL.t = var->t;
	SL.n = n;

	triang(&SL);
	retrossubs(&SL, var->delta);
}

