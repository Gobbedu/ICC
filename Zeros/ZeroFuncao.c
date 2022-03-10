#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"


double bisseccao (Polinomio p, double a, double b, double eps,
	       int *it, double *raiz)
{
	double EAr, m_old, m_new = 0;
	double Pxa, Dxa, Pxm, Dxm;
	int iter = 0;

	// printf("iter \t m_old \t\t m_new \t\t a \t\t\t b \t\t Pxa \t\t Pxb\n");
	// printf("---- \t ----- \t\t ----- \t\t ----- \t\t ----- \t\t ----- \t\t -----\n");

	do {
		m_old = m_new;
		m_new = (a + b)/2;

		calcPolinomio_rapido(p, a, &Pxa, &Dxa);
		calcPolinomio_rapido(p, m_new, &Pxm, &Dxm);

		if( Pxa * Pxm < 0 )
			b = m_new;
		else if (Pxa * Pxm > 0)
			a = m_new;
		else{
			*raiz = m_new;
			*it = iter;
			return 0;
		}

		EAr = fabs(m_new - m_old);
		iter++;

	} while ( iter < MAXIT && EAr > eps);

	*it = iter;
	*raiz = m_new;

	// o valor a retornar eh o erro da raiz calculado
	return EAr;
		
}


double newtonRaphson(Polinomio p, double x0, double eps,
		   int *it, double *raiz)
{
	double EAr, Px, Dx, x_old, x_new = x0;
	double fnew, fold;
	int crit1, crit2, crit3;
	int iter = 0;

	// x_new = phi(x_old);

	do{
		x_old = x_new;

		calcPolinomio_rapido(p, x_old, &Px, &Dx);

		// fk1 = k0 - fk0/f'k0 
		x_new = x_old - (Px/Dx); 

		calcPolinomio_rapido(p, x_new, &fnew, &Dx);
		// calcPolinomio_rapido(p, x_old, &fold, &Dx);

		crit1 = (fabs(x_old - x_new) > eps);
		crit2 = (fabs(fnew) > eps);
		crit3 = (iter < MAXIT);
		iter++;
	} while( crit1 && crit2 && crit3 );

	*it = iter;
	*raiz = x_new;
	 
	 // returns error of root
	 return fabs((x_new - x_old)/ x_new) * 100;

}


double secante (Polinomio p, double x0, double x1, double eps,
	     int *it, double *raiz)
{
	double EAr, Px, Dx, x_oold, x_old, x_new;
	double fnew, fold, foold;
	int crit1, crit2, crit3, iter = 0;
	x_old = x0;
	x_new = x1;


	do{
		x_oold = x_old;
		x_old = x_new;

		calcPolinomio_rapido(p, x_old, &fold, &Dx);
		calcPolinomio_rapido(p, x_oold, &foold, &Dx);

		// x_new = phi(x_old);
		x_new = x_old - ((fold*(x_old - x_oold))/(fold - foold));

		iter++;

		EAr = fabs((x_new - x_old)/ x_new) * 100;

		crit1 = (fabs(x_old - x_new) > EPS);
		crit2 = (fabs(fnew - fold) > DBL_EPSILON);
		crit3 = (iter < MAXIT);
		
		// printf("crit1: %f < %f => %d\n", fabs(x_old - x_new), EPS, crit1);
		// printf("crit2: |%f - %f| < %f => %d\n",fnew, fold, EPS, crit2);
		// printf("EAr: %1.9ef  EPS: %1.9ef\n", EAr, eps);
		
	} while( crit1 && crit2 && crit3 );

	*it = iter;
	*raiz = x_new;
	 
	 // returns error of root
	 return EAr;
}


void calcPolinomio_rapido(Polinomio p, double x, double *Px, double *Dx)
{
	double b = 0;
	double c = 0;
	
	for(int i = p.grau; i > 0; i--)
	{
		b = b*x + p.p[i];
		c = c*x + b;
	}
	b = b*x + p.p[0];

	// printf("resultado calcPol: Px = %f\n", b);
	// printf("resultado calcPol:Dx = %f\n", c);
	*Px = b;
	*Dx = c;
}


void calcPolinomio_lento(Polinomio p, double x, double *Px, double *Dx)
{
	// TODO libmatheval

	for(int i = p.grau ; i >= 0; i--)
	{
		*Px += p.p[i]*pow(x, i);
	}
	
	*Dx = 0;
}
