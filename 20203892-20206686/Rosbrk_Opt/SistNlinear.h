/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#ifndef __SIST_N_LINEAR__
#define __SIST_N_LINEAR__

#include "utils.h"

#define STR_BUFFER 150000

typedef struct {
    int n;
    double *chute;
    double eps;
    int iteracao;
    char funcao[STR_BUFFER];

    void *f;        // funcao original
    void ***Hf;     // Hessiana de funcoes
    void **Gf;      // Gradiente de funcoes
    char **names;   // nomes variaveis (para matheval)
} SistNl_t;

typedef struct {
    double *x0;
    double *x1;
    double *delta;

    // SistLinear_t *sl;
    int *t;

    double **He;    // Hessiana exata (com valores)
    double *Ge;     // Gradiente exata (com valores)
} SnlVar_t;

typedef struct {
  unsigned int n;   // dimensão do SL
  double **A;       // matriz dos coeficientes do SL (vetor de ponteiros para posições de M)
  double *b;        // termos independentes do SL
  int *t;           // trocas efetuadas em LU
} SistLinear_t;


SistNl_t *alocaSistNl(unsigned int n);
SnlVar_t *alocaSnlVar(double *chute, int n);
SistNl_t *lerSistNL(void);


void genGradiente(SistNl_t *snl);
void genHessiana(SistNl_t *snl);
void genNames(SistNl_t *snl);

void substituteX(SistNl_t *snl, SnlVar_t *nt);
void calcHessiana(double **HessExata, int n, double *x0); 
void calcGradiente(double *GradExato, int n, double *x0);
void calcDelta(SnlVar_t *var, int tam);

int Parada(double *delta, double eps, int n);

void snl2sl(SistNl_t *snl, SnlVar_t *nt);
void printCol(double* pto, int i, int max, FILE *saida);
void varinfo(SnlVar_t var, SistNl_t snl);
void snlinfo(SistNl_t *snl);

void liberaSnlVar(SnlVar_t *var, int n);
void liberaMatheval(SistNl_t *snl);
void liberaSistNl(SistNl_t *SL);

#endif
