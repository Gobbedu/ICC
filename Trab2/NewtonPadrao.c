/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "NewtonPadrao.h"
// #include <likwid.h>
// // ELIMINACAO GAUSS ou NEWTON PADRAO
void NewtonPadrao(SistNl_t *snl, double* resposta, Tempo_t *t, int *nIter)
{    
	// inicializa tudo que precisa
	// double *resposta = malloc(sizeof(double) * snl->iteracao);
	double tTotal, tauxder, tGrad, tHess, tauxSL;

	int itr = 0;


	SnlVar_t *np = alocaSnlVar(snl->chute, snl->n);

	// LIKWID_MARKER_INIT;

	// LIKWID_MARKER_START("marker-METODO");
	// --------LOOP PRINCIPAL-------- //
	for(int i = 0; i < snl->iteracao; i++)
	{
		tTotal = timestamp();

		// tauxder = timestamp();
		// LIKWID_MARKER_START("marker-GRADIENTE");
		tGrad = timestamp();
		calcGradiente(snl, np);					// calcula J[X]
		tGrad = timestamp() - tGrad;
		// LIKWID_MARKER_STOP("marker-GRADIENTE");			

		// LIKWID_MARKER_START("marker-HESSIANA");
		tHess = timestamp();
		calcHessiana(snl, np);					// calcula H[X]
		tHess = timestamp() - tHess;
		// LIKWID_MARKER_STOP("marker-HESSIANA");			
		// tauxder = timestamp() - tauxder;

		snl2sl(snl, np);                        // copia dados de snl em sl

		#ifndef ROSENBROCK
			resposta[i] = evaluator_evaluate(snl->f, snl->n, snl->names, np->x0);
		#else
			resposta[i] = rosenbrock(np->x0, snl->n);
		#endif
		
		// LIKWID_MARKER_START("marker-SIST-LINEAR");
		tauxSL = timestamp();
		eliminacaoGauss(np->sl, np->delta);     // calcula H[X]*delta = - J[X]  // A*x = -b
		tauxSL = timestamp() - tauxSL;
		// LIKWID_MARKER_STOP("marker-SIST-LINEAR");			

		// varinfo(*np, *snl);
		calcDelta(np, snl->n);                     // X[i+1] = X[i] + delta[i]

		// Devolve f(X), ponto critico estimado

		tTotal = timestamp() - tTotal;

		// t->totalMetodo += tTotal;
		// t->Gradiente += tGrad;
		// t->Hessiana += tHess;
		// t->totalSL += tauxSL;
		// t->derivadas += tauxder;
		itr++;
		
			if(Parada(snl, np->delta))
			break;
	}
	// LIKWID_MARKER_STOP("marker-METODO");			

	// LIKWID_MARKER_CLOSE;

	*nIter = itr;

	liberaSnlVar(np, snl->n);
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

