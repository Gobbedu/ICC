/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/

#include "SistNlinear.h"

SistNl_t *lerSistNL(void)
{
    unsigned int n;

    SistNl_t *SnL = NULL;
  
    if (scanf("%d",&n) != EOF) 
    {  
        SnL = alocaSistNl(n);
        if (!SnL) return NULL;
        
        assert(scanf("%s", SnL->funcao) > 0);

        for(int i = 0; i < n; i++)
            assert(scanf("%lf", &(SnL->chute[i])) > 0);

        assert(scanf("%lf", &(SnL->eps)) > 0);
        assert(scanf("%i", &(SnL->iteracao)) > 0);

    }
  
  return SnL;
}

int Parada(SistNl_t *snl, double *delta)
{
    // normal -> || X || = max{ |xi|, 1 ≤ i ≤ n}
    double maxD = -INFINITY;
    for(int i = 0; i < snl->n; i++){
        maxD = (fabs(delta[i]) > maxD) ? fabs(delta[i]) : maxD;
    }

    // se ||delta|| < eps PARE // se ||Ji(X)|| < eps PARE
    return (maxD < snl->eps);
}

void calcDelta(SnlVar_t *var, int n){
    // allow for loop unrolling
    for(int i = 0; i < n; i++)
    {
        var->x1[i] = var->x0[i] + var->delta[i];
        var->x0[i] = var->x1[i];
    }

}

SistNl_t *alocaSistNl(unsigned int n){
    SistNl_t *SnL = (SistNl_t *) malloc(sizeof(SistNl_t));

    if (SnL){
        SnL->n = n;

        // vetor de chutes
        SnL->chute = malloc(sizeof(double)*n);

        // snl->names[i] = malloc(4 * sizeof(char)); // maximo 999 variaveis
        SnL->names = malloc(SnL->n * sizeof(char *));
        for(int i = 0; i < SnL->n; i++)
            SnL->names[i] = malloc(5 * sizeof(char));
    }
  return SnL;
}

void liberaSistNl(SistNl_t *snl) {
    for(int i = 0; i < snl->n; i++)
    {
        free(snl->names[i]);
    }

    free(snl->names);
    free(snl->chute);
    free(snl);
}

SnlVar_t *alocaSnlVar(double *chute, int n)
{
    SnlVar_t *var = malloc(sizeof(SnlVar_t));
    var->delta = malloc(n * sizeof(double));
    var->x0 = malloc(n * sizeof(double));
    var->x1 = malloc(n * sizeof(double));
    var->sl = alocaSistLinear(n);

    for(int i = 0; i < n; i++)
    {
        var->x0[i] =    0;
        var->x1[i] =    0;
        var->delta[i] = 1;
    }

    // aloca He e libera var caso erro
    var->He = (double **) malloc(sizeof(double*)*n);
    if (!(var->He)) {
        liberaSistLinear(var->sl);
        free(var->x0);
        free(var->x1);
        free(var->delta);
        free(var);
        fprintf(stderr, "erro ao alocar Hessiana exata em variaveis do sist_n_linear\n");
        return NULL;
    }
    for(int i = 0; i < n; i++){
        var->He[i] = malloc(sizeof(double)*n);
    }
        
    var->Ge = (double *) malloc(sizeof(double)*n);
    if (!(var->Ge)) {
        liberaSistLinear(var->sl);
        free(var->x0);
        free(var->x1);
        free(var->delta);
        free(var->He);
        free(var);
        fprintf(stderr, "erro ao alocar Gradiente exata em variaveis do sist_n_linear\n");
        return NULL;
    }

    // chute inicial
    for(int i = 0; i < n; i++) 
        var->x0[i] = chute[i]; 

    return var;
}

void liberaSnlVar(SnlVar_t *var, int n)
{
    liberaSistLinear(var->sl);
    free(var->x0);
    free(var->x1);
    free(var->delta);

    for(int i = 0; i < n; i++)
    {
        free(var->He[i]);
    }
    free(var->He);
    free(var->Ge);

    free(var);
}

void calcHessiana(SistNl_t *snl, SnlVar_t *var)
{
    for(int i = 0; i < snl->n ; i++)
        for(int j = 0; j < snl->n; j++)
            var->He[i][j] = rosenbrock_dxdy(i, j, var->x0, snl->n);
}

void calcGradiente(SistNl_t *snl, SnlVar_t *var)
{
    for(int i = 0; i < snl->n ; i++)
        var->Ge[i] = rosenbrock_dx(i, var->x0, snl->n);
}
/*
void substituteX(SistNl_t *snl, SnlVar_t *nt)
{
    for(int i = 0; i < snl->n ; i++)
    {
        for(int j = 0; j < snl->n; j++)
            nt->He[i][j] = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, nt->x0);// H[X]
        nt->Ge[i] = evaluator_evaluate(snl->Gf[i], snl->n, snl->names, nt->x0);// J[X]
    }
}
*/

void snl2sl(SistNl_t *snl, SnlVar_t *nt)
{
    for(int i = 0; i < snl->n; ++i) {
        for(int j = 0; j < snl->n; ++j)
        {
            nt->sl->A[i][j] = nt->He[i][j];
        }
        nt->sl->b[i] = - nt->Ge[i];
    }
}
/*
void snl2sl(SistNl_t *snl, SnlVar_t *nt)
{
    nt->sl->A = nt->He;
    nt->sl->b = nt->Ge;
}
*/

void printCol(double* pto, int i, int max, FILE *saida)
{
    if(i < max){
        if (isnan(pto[i]) || isinf(pto[i]))
            fprintf(saida, "%1.14e\t\t\t| ", pto[i]);
        else
            fprintf(saida, "%1.14e\t| ", pto[i]);
    }
    else
        fprintf(saida, "\t\t\t| ");
}
