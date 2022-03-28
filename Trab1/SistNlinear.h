#ifndef __SIST_N_LINEAR__
#define __SIST_N_LINEAR__

#ifndef __SIST_LINEAR__
#include "SistLinear.h"
#endif

#define HESS_STEP 1
#define STR_BUFFER 1000

typedef struct {
        SistLinear_t *sl;
        double *x0;
        double *x1;
        double *delta;
} SnlVar_t;

typedef struct {
    int n;
    double *chute;
    double eps;
    int iteracao;
    char funcao[STR_BUFFER];

    char **names;
    
    void *f;        // funcao original

    void ***Hf;     // Hessiana de funcoes
    double **He;    // Hessiana exata (com valores)

    void **Bf;      // Jacobiana de funcoes
    double *Be;     // Jacobiana exata (com valores)

} SistNl_t;


SistNl_t *alocaSistNl(unsigned int n);
SnlVar_t *genSnlVar(double *chute, int n);
SistNl_t *lerSistNL(void);

double NewtonPadrao(SistNl_t *snl, SnlVar_t *np);
double NewtonModificado(SistNl_t *snl, SnlVar_t *nm);
double NewtonInexato(SistNl_t *snl, SnlVar_t *ni);

double *genValues(int n, double init);
void genSistNaoLinear(SistNl_t *snl);
void genJacobiana(SistNl_t *snl);
void liberaSistNl(SistNl_t *SL);
void genHessiana(SistNl_t *snl);
void genNames(SistNl_t *snl);

void printCol(double pto, SistNl_t *snl, SnlVar_t *nt);

void calcDelta(SnlVar_t *var, int n);
void substituteX(SistNl_t *snl, double *X);
double minDelta(double *delta, int n);
void snlinfo(SistNl_t *snl);
void varinfo(SnlVar_t var, SistNl_t snl);

void liberaMatheval(SistNl_t *snl);
void snl2sl(SistNl_t *snl, SistLinear_t *sl);
void liberaSnlVar(SnlVar_t *var);


#endif
