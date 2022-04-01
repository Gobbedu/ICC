#ifndef __SIST_N_LINEAR__
#define __SIST_N_LINEAR__

#ifndef __SIST_LINEAR__
#include "SistLinear.h"
#endif

#define HESS_STEP 1
#define STR_BUFFER 1000

typedef struct {
    int n;
    double *chute;
    double eps;
    int iteracao;
    char funcao[STR_BUFFER];

    void *f;        // funcao original
    void ***Hf;     // Hessiana de funcoes
    void **Bf;      // Jacobiana de funcoes
    char **names;   // nomes variaveis (para matheval)
} SistNl_t;

typedef struct {
    double *x0;
    double *x1;
    double *delta;

    SistLinear_t *sl;

    double **He;    // Hessiana exata (com valores)
    double *Je;     // Jacobiana exata (com valores)
} SnlVar_t;


SistNl_t *alocaSistNl(unsigned int n);
SnlVar_t *alocaSnlVar(double *chute, int n);
SistNl_t *lerSistNL(void);

double NewtonModificado(SistNl_t *snl, SnlVar_t *nm, int i);
double NewtonInexato(SistNl_t *snl, SnlVar_t *ni);

double *genValues(int n, double init);
void genSistNaoLinear(SistNl_t *snl);
void genJacobiana(SistNl_t *snl);
void genHessiana(SistNl_t *snl);
void genNames(SistNl_t *snl);

void substituteX(SistNl_t *snl, SnlVar_t *nt);
void calcDelta(SistNl_t *snl, SnlVar_t *var);
void printCol(double* pto, int i, int max);
void snl2sl(SistNl_t *snl, SnlVar_t *nt);
void varinfo(SnlVar_t var, SistNl_t snl);
int Parada(SistNl_t *snl, double *delta);
void snlinfo(SistNl_t *snl);

void liberaSnlVar(SnlVar_t *var, int n);
void liberaMatheval(SistNl_t *snl);
void liberaSistNl(SistNl_t *SL);


#endif
