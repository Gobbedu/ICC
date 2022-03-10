#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "utils.h"
#include "ZeroFuncao.h"


void root_eval(Polinomio pol, double a, double b, double x0, double x1, double eps){
  int it;
  double raiz, erro, t1,t2;
  
  it = raiz = erro = 0;
  t1 = timestamp();
  erro = bisseccao(pol, a, b, eps, &it, &raiz);
  t2 = timestamp();
  printf("bisseccao  \t %d \t \t %1.9ef \t %1.9ef \t%1.5ef\n", it, raiz, erro, t2-t1);


  it = raiz = erro = 0;
  t1 = timestamp();
  erro = newtonRaphson(pol, x0, eps, &it, &raiz);
  t2 = timestamp();
  printf("newton  \t %d \t \t %1.9ef \t %1.9ef \t%1.5ef\n", it, raiz, erro, t2-t1);


  it = raiz = erro = 0;
  t1 = timestamp();
  erro = secante(pol, x0, x1, eps, &it, &raiz);
  t2 = timestamp();
  printf("secante  \t %d \t \t %1.9ef \t %1.9ef \t%1.5ef\n\n", it, raiz, erro, t2-t1);

}


int main ()
{
  Polinomio pol;
  double eps = EPS;


  printf("\nmetodo \t\t iteracoes \t raiz \t\t\t erro \t\t\ttime \n");
  printf("--------- \t --------- \t ----------------\t ---------------- \t------------\n");

  // 2x⁴ + 4x³ + 3x² − 10x − 15 = 0,    ξ ∈ [0, 3] #root -> 1.4928
  double p1[5] = {-15, -10, 3, 4, 2};
  pol.grau = 4;
  pol.p = p1;
  root_eval(pol, 0, 3, 2, 2.5, eps);


  // x⁵ − 2x⁴ − 9x³ + 22x² + 4x − 24 = 0,  ξ ∈ [0, 5] #root -> 2.0 \ -3 \ -1
  double p2[6] = {-24, 4, 22, -9, -2, 1};
  pol.grau = 5;
  pol.p = p2;
  root_eval(pol, 1, 2, 0, 1, eps);


  // 3x⁴ + 2x³ + 4x² + 25x − 30 = 0,      ξ ∈ [0, 1]  #root -> 0.91803
  double p3[5] = {-30, 25, 4, 2, 3};
  pol.grau = 5;
  pol.p = p3;
  root_eval(pol, 0, 5, 0, 0.5, eps);


  // 2x³ − 5x² − x + 3 = 0,           ξ ∈ [0, 1]    #root -> 0.80465
  double p4[4] = {3, -1, -5, 2};
  pol.grau = 3;
  pol.p = p4;
  root_eval(pol, 0, 1, 0.3, 0.5, eps);

  printf("\n");


  void *f;
  f = evaluator_create("x");
  assert(f);
  char *names[] = { "x" };
  double values[] = { 0.1 };
  printf ("f(0.1) = %g\n", evaluator_evaluate (f, 1, names,
                                                    values));

  return 0;
}

