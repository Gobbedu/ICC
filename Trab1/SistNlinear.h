/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#ifndef __SIST_N_LINEAR__
#define __SIST_N_LINEAR__

#include "utils.h"
#include "SistLinear.h"

#define STR_BUFFER 1000

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

    SistLinear_t *sl;

    double **He;    // Hessiana exata (com valores)
    double *Ge;     // Gradiente exata (com valores)
} SnlVar_t;


SistNl_t *alocaSistNl(unsigned int n);
SnlVar_t *alocaSnlVar(double *chute, int n);
SistNl_t *lerSistNL(void);


void genGradiente(SistNl_t *snl);
void genHessiana(SistNl_t *snl);
void genNames(SistNl_t *snl);

void substituteX(SistNl_t *snl, SnlVar_t *nt);
void calcHessiana(SistNl_t *snl, SnlVar_t *var); 
void calcGradiente(SistNl_t *snl, SnlVar_t *var);
void calcDelta(SnlVar_t *var, int tam);
int Parada(SistNl_t *snl, double *delta);


void snl2sl(SistNl_t *snl, SnlVar_t *nt);
void printCol(double* pto, int i, int max);
void varinfo(SnlVar_t var, SistNl_t snl);
void snlinfo(SistNl_t *snl);

void liberaSnlVar(SnlVar_t *var, int n);
void liberaMatheval(SistNl_t *snl);
void liberaSistNl(SistNl_t *SL);

#endif
