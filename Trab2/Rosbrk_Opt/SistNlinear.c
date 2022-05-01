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

int Parada(double *delta, double eps, int n)
{
    // normal -> || X || = max{ |xi|, 1 ≤ i ≤ n}
    double maxD = -INFINITY;
    for(int i = 0; i < n && maxD < eps; i++){
        maxD = (fabs(delta[i]) > maxD) ? fabs(delta[i]) : maxD;
    }

    // se ||delta|| < eps PARE // se ||G(X)|| < eps PARE
    return (maxD < eps);
}

void calcDelta(SnlVar_t *var, int n){
    // allow for loop unrolling
    for(int i = 0; i < n - (n%8); i+=8){
        var->x1[i] = var->x0[i] + var->delta[i];
        var->x1[i+1] = var->x0[i+1] + var->delta[i+1];
        var->x1[i+2] = var->x0[i+2] + var->delta[i+2];
        var->x1[i+3] = var->x0[i+3] + var->delta[i+3];
        var->x1[i+4] = var->x0[i+4] + var->delta[i+4];
        var->x1[i+5] = var->x0[i+5] + var->delta[i+5];
        var->x1[i+6] = var->x0[i+6] + var->delta[i+6];
        var->x1[i+7] = var->x0[i+7] + var->delta[i+7];
    }
    for(int i = n-(n%8); i < n; i++)
        var->x1[i] = var->x0[i] + var->delta[i];

    for(int i = 0; i < n-(n%8); i+=8){
        var->x0[i] = var->x1[i];
        var->x0[i+1] = var->x1[i+1];
        var->x0[i+2] = var->x1[i+2];
        var->x0[i+3] = var->x1[i+3];
        var->x0[i+4] = var->x1[i+4];
        var->x0[i+5] = var->x1[i+5];
        var->x0[i+6] = var->x1[i+6];
        var->x0[i+7] = var->x1[i+7];
    }
    for(int i = n - (n%8); i < n; i++)
        var->x0[i] = var->x1[i];

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
    // var->sl = alocaSistLinear(n);
    var->t = malloc(n * sizeof(int));

    memset(var->x0, 0, n*sizeof(double));
    memset(var->x1, 0, n*sizeof(double));
    memset(var->delta, 1, n*sizeof(double));

    for(int i = 0; i < n; i++)
        var->t[i] = i;

    // aloca He e libera var caso erro
    var->He = (double **) malloc(sizeof(double*)*n);
    if (!(var->He)) {
        // liberaSistLinear(var->sl);
        free(var->delta);
        free(var->x0);
        free(var->x1);
        free(var->t);
        free(var);
        fprintf(stderr, "erro ao alocar Hessiana exata em variaveis do sist_n_linear\n");
        return NULL;
    }

    int padding;

    //regra de ouro
    if( (n & (n - 1)) == 0) padding = 1;
    else padding = 0;

    for(int i = 0; i < n; i++){
        var->He[i] = malloc(sizeof(double)*n + padding);
    }
        
    var->Ge = (double *) malloc(sizeof(double)*n);
    if (!(var->Ge)) {
        // liberaSistLinear(var->sl);
        free(var->delta);
        free(var->He);
        free(var->x0);
        free(var->x1);
        free(var->t);
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
    // liberaSistLinear(var->sl);
    // free(var->t);
    free(var->x0);
    free(var->x1);
    free(var->delta);

    for(int i = 0; i < n; i++){
        free(var->He[i]);
    }

    free(var->t);

    free(var->He);
    free(var->Ge);
    free(var);
}

void calcHessiana(double **He, int n, double *x0)
{
    for(int i = 0; i < n ; i++){
        for(int j = 0; j < n-(n%8); j+=8){
            He[i][j] = rosenbrock_dxdy(i, j, x0, n);
            He[i][j+1] = rosenbrock_dxdy(i, j+1, x0, n);
            He[i][j+2] = rosenbrock_dxdy(i, j+2, x0, n);
            He[i][j+3] = rosenbrock_dxdy(i, j+3, x0, n);
            He[i][j+4] = rosenbrock_dxdy(i, j+4, x0, n);
            He[i][j+5] = rosenbrock_dxdy(i, j+5, x0, n);
            He[i][j+6] = rosenbrock_dxdy(i, j+6, x0, n);
            He[i][j+7] = rosenbrock_dxdy(i, j+7, x0, n);
        }
        for(int j = n-(n%8); j < n; j++)
            He[i][j] = rosenbrock_dxdy(i, j, x0, n);
    }
}

void calcGradiente(double *Ge, int n, double *x0)
{
    for(int i = 0; i < n-(n%8); i+=8){
        Ge[i] = rosenbrock_dx(i, x0, n);
        Ge[i+1] = rosenbrock_dx(i+1, x0, n);
        Ge[i+2] = rosenbrock_dx(i+2, x0, n);
        Ge[i+3] = rosenbrock_dx(i+3, x0, n);
        Ge[i+4] = rosenbrock_dx(i+4, x0, n);
        Ge[i+5] = rosenbrock_dx(i+5, x0, n);
        Ge[i+6] = rosenbrock_dx(i+6, x0, n);
        Ge[i+7] = rosenbrock_dx(i+7, x0, n);
    }
    for(int i = 0; i < n ; i++)
        Ge[i] = rosenbrock_dx(i, x0, n);
}

/*
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
