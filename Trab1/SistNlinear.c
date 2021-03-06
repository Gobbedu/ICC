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
        
        scanf("%s", SnL->funcao);

        for(int i = 0; i < n; i++)
            scanf("%lf", &(SnL->chute[i]));

        scanf("%lf", &(SnL->eps));
        scanf("%i", &(SnL->iteracao));

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
        var->x1[i] = var->x0[i] + var->delta[i];

    for(int i = 0; i < n; i++)
        var->x0[i] = var->x1[i];

}

void genNames(SistNl_t *snl){  
    for(int i = 0; i < snl->n; i++){
        sprintf(snl->names[i], "x%i", i+1);
    }
}

void genGradiente(SistNl_t *snl){
    for(int i = 0; i < snl->n; i++){
        snl->Gf[i] = evaluator_derivative(snl->f, snl->names[i]);
    }
}

void genHessiana(SistNl_t *snl){
    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            // deriva Jacob[i] n vezes em x1..xn
            snl->Hf[i][j] = evaluator_derivative(snl->Gf[i], snl->names[j]);
        }
    }
}

SistNl_t *alocaSistNl(unsigned int n){
    SistNl_t *SnL = (SistNl_t *) malloc(sizeof(SistNl_t));

    if (SnL){
        SnL->n = n;

        // vetor de chutes
        SnL->chute = malloc(sizeof(double)*n);

        // snl->names[i] = malloc(4 * sizeof(char)); // maximo 9999 variaveis
        SnL->names = malloc(SnL->n * sizeof(char *));
        for(int i = 0; i < SnL->n; i++)
            SnL->names[i] = malloc(5 * sizeof(char));

        SnL->Hf = (void ***) malloc(sizeof(void **)*n);
        if (!(SnL->Hf)) {
            free(SnL);
            return NULL;
        }
        for (int i=0; i < n; ++i)
            SnL->Hf[i] = (void **) malloc(sizeof(void *)*n);


        SnL->Gf = (void **) malloc(sizeof(void *)*n);
        if(!(SnL->Gf)){
            free(SnL->Hf);
            free(SnL);
        }

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

void liberaMatheval(SistNl_t *snl)
{
    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            evaluator_destroy(snl->Hf[i][j]);
        }
        evaluator_destroy(snl->Gf[i]);
        free(snl->Hf[i]);
    }
    evaluator_destroy(snl->f);
    free(snl->Gf);
    free(snl->Hf);
}

SnlVar_t *alocaSnlVar(double *chute, int n)
{
    // CALLOC 
    SnlVar_t *var = malloc(sizeof(SnlVar_t));
    var->delta = malloc(n * sizeof(double));
    var->x0 = malloc(n * sizeof(double));
    var->x1 = malloc(n * sizeof(double));
    var->sl = alocaSistLinear(n);

    // VETOR UNICO DENTRO DO FOR
    for(int i = 0; i < n; i++)
    {
        var->x0[i] =    0;
        var->x1[i] =    0;
        var->delta[i] = 1;
    }
    // TRAB2: VETORZAO ->
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

// unico for p/ calc hessiana
void calcHessiana(SistNl_t *snl, SnlVar_t *var)
{
    for(int i = 0; i < snl->n ; i++)
        for(int j = 0; j < snl->n; j++)
            var->He[i][j] = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, var->x0);
}

void calcGradiente(SistNl_t *snl, SnlVar_t *var)
{
    for(int i = 0; i < snl->n ; i++)
        var->Ge[i] = evaluator_evaluate(snl->Gf[i], snl->n, snl->names, var->x0);
}

// ruim, nn usar
void substituteX(SistNl_t *snl, SnlVar_t *nt)
{
    for(int i = 0; i < snl->n ; i++)
    {
        for(int j = 0; j < snl->n; j++)
            nt->He[i][j] = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, nt->x0);// H[X]
        nt->Ge[i] = evaluator_evaluate(snl->Gf[i], snl->n, snl->names, nt->x0);// J[X]
    }
}

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


void snlinfo(SistNl_t *S){
    printf("\n------------------------------------SNL INFO------------------------------------\n");
    printf("#  n: %i\n", S->n);
    printf("#  funcao: %s\n", S->funcao);

    printf("#  chute:");
    for(int i = 0; i < S->n; i++)
        printf(" %g ", S->chute[i]);
    printf("\n");

    printf("#  eps: %g\n", S->eps);
    printf("#  iteracao: %i\n", S->iteracao);

    printf("#  names: ");
    for(int i = 0; i < S->n ; i++)
        printf(" %s ", S->names[i]);
    printf("\n");

    printf("\nVETOR CHUTE's:\n#");
    for(int i = 0; i < S->n; i++)
        printf("  %f, ", S->chute[i]);
    printf("\nEVALUATOR\n#  f(CHUTE) = %g\n",evaluator_evaluate(S->f, S->n, S->names, S->chute));


    printf("--------------------------------------END SNL INFO-----------------------------------------\n");
}

void varinfo(SnlVar_t nt, SistNl_t snl)
{
    char *func;
    printf("\n------------------------------------INFO VARIAVEIS------------------------------------\n");
    printf("x0: "); for(int i = 0; i < snl.n; i++) printf("%1.3e ", nt.x0[i]); printf("\n");
    printf("x1: "); for(int i = 0; i < snl.n; i++) printf("%1.3e ", nt.x1[i]); printf("\n");
    printf("delta: "); for(int i = 0; i < snl.n; i++) printf("%1.3e ", nt.delta[i]); printf("\n");
    printf("evaluator[x0] = %1.14e\n", evaluator_evaluate(snl.f, snl.n, snl.names, nt.x0));
    
        

    printf("-------------------------------------------------------------------------------\n");
}

void printCol(double* pto, int i, int max, int argc, FILE *saida)
{
    if(argc == 3){
        if(i < max){
            if (isnan(pto[i]) || isinf(pto[i]))
                fprintf(saida, "%1.14e\t\t\t| ", pto[i]);
            else
                fprintf(saida, "%1.14e\t| ", pto[i]);
        }
        else
            fprintf(saida, "\t\t\t| ");
    }
    else{
        if(i < max){
            if (isnan(pto[i]) || isinf(pto[i]))
                printf("%1.14e\t\t\t| ", pto[i]);
            else
                printf("%1.14e\t| ", pto[i]);
        }
        else
            printf("\t\t\t| ");
    }

}
