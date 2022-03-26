#ifndef __SIST_N_LINEAR__
#define __SIST_N_LINEAR__

#define HESS_STEP 1
#define STR_BUFFER 1000


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
void liberaSistNl(SistNl_t *SL) ;
SistNl_t *lerSistNL(void);

void genNames(SistNl_t *snl);
double *genValues(int n, double init);
void genJacobiana(SistNl_t *snl);
void genHessiana(SistNl_t *snl);
void genSistNaoLinear(SistNl_t *snl);

void substituteX(SistNl_t *snl, double *X);
void calcDelta(double *new_values, double *old_values, double *delta, int n);
double minDelta(double *delta);
void snlinfo(SistNl_t *snl);


#endif
