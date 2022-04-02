
#include "SistNlinear.h"

/*Gauss seidel do Vods, n funciona 100%*/
void gauss_seidel(SistLinear_t *SL, double *X)
{
    int i,j,q,d;
  unsigned int n = SL->n;
  double *r =  malloc((SL->n) * sizeof (double)) ;
  double temp,sum,erroMaximo,erroCalculado;
  double **A = SL->A;
  double *b = SL->b;
  double *x = malloc((SL->n) * sizeof(double));
  //X = malloc((SL->n) * sizeof(double));

  // A[n][n] = Matriz principal (e->f)
  // b[n] = vetor_independente (e->termos_independentes)

  for(i=0;i<n;i++){
      r[i] = 0;
  }
  
  q = 0;
  do{
      erroCalculado = 0;
      q++;
      for(i=0;i<n;i++){
          sum = 0;
          for(j=0;j<n;j++){
              if(i != j){
                  sum = sum + (A[i][j] * r[j]);
              }
            }
          temp = (-1.0 / A[i][i]) * sum + b[i] / A[i][i];
          erroMaximo = fabs(temp - r[i]);
          r[i] = temp;
          if(erroMaximo > erroCalculado)
              erroCalculado = erroMaximo;
        }
  }while(erroCalculado >= 1e-6  && q<=50);

  for(i=0;i<n;i++){ //copiar dados calculados para *x
    x[i] = r[i];
  }

  //X= x; //copiar dados para estrutura

  free(r);
  free(x);
  //free(X);


}
double NewtonInexato(SistNl_t *snl, SnlVar_t *ni)
{
    /* NEWTON INEXATO

    // GAUSS SEIDEL (TODO)

    */

    //PODE DELETAR AKI
    substituteX(snl, ni);               // calcula H[X] e J[X]
    // printCol(evaluator_evaluate(snl->f, snl->n, snl->names, ni), snl, ni);
    snl2sl(snl, ni);                    // copia dados de snl em sl
    gauss_seidel(ni->sl,ni->delta);         // calcula H[X]*delta = - J[X]  // A*x = -b
    calcDelta(snl, ni);                  // X[i+1] = X[i] + delta[i]
    return evaluator_evaluate(snl->f, snl->n, snl->names, ni->x1);  

}


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
    // normal -> || J(X) || = max{ |Ji(X)|, 1 ≤ i ≤ n}
    // double maxF = -INFINITY;
    // double Ji;
    // for(int i = 0; i < snl->n; i++)
    // {
    //     Ji = evaluator_evaluate(snl->Bf[i], snl->n, snl->names, nt->x0);
    //     maxF = (fabs(Ji) > maxF) ? fabs(Ji) : maxF;
    // }

    // normal -> || X || = max{ |xi|, 1 ≤ i ≤ n}
    double maxD = -INFINITY;
    for(int i = 0; i < snl->n; i++){
        maxD = (fabs(delta[i]) > maxD) ? fabs(delta[i]) : maxD;
    }

    // se ||delta|| < eps PARE // se ||Ji(X)|| < eps PARE
    return (maxD < snl->eps);
}

void calcDelta(SistNl_t *snl, SnlVar_t *var){
    // if(PARADA)
    for(int i = 0; i < snl->n; i++){
        var->x1[i] = var->x0[i] + var->delta[i];
        var->x0[i] = var->x1[i];
    }
}

void genNames(SistNl_t *snl){  
    for(int i = 0; i < snl->n; i++){
        sprintf(snl->names[i], "x%i\0", i+1);
    }
}

void genJacobiana(SistNl_t *snl){
    for(int i = 0; i < snl->n; i++){
        snl->Bf[i] = evaluator_derivative(snl->f, snl->names[i]);
    }
}

void genHessiana(SistNl_t *snl){
    for(int i = 0; i < snl->n; i++){
        for(int j = 0; j < snl->n; j++)
        {
            // deriva Jacob[i] n vezes em x1..xn
            snl->Hf[i][j] = evaluator_derivative(snl->Bf[i], snl->names[j]);
        }
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

        SnL->Hf = (void ***) malloc(sizeof(void **)*n);
        if (!(SnL->Hf)) {
            free(SnL);
            return NULL;
        }
        for (int i=0; i < n; ++i)
            SnL->Hf[i] = (void **) malloc(sizeof(void *)*n);


        SnL->Bf = (void **) malloc(sizeof(void *)*n);
        if(!(SnL->Bf)){
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
        evaluator_destroy(snl->Bf[i]);
        free(snl->Hf[i]);
    }
    evaluator_destroy(snl->f);
    free(snl->Bf);
    free(snl->Hf);
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
        
    var->Je = (double *) malloc(sizeof(double)*n);
    if (!(var->Je)) {
        liberaSistLinear(var->sl);
        free(var->x0);
        free(var->x1);
        free(var->delta);
        free(var->He);
        free(var);
        fprintf(stderr, "erro ao alocar Jacobiana exata em variaveis do sist_n_linear\n");
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
    free(var->Je);

    free(var);
}

void calcHessiana(SistNl_t *snl, SnlVar_t *var)
{
    for(int i = 0; i < snl->n ; i++)
        for(int j = 0; j < snl->n; j++)
            var->He[i][j] = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, var->x0);
}

calcJacobiana(SistNl_t *snl, SnlVar_t *var)
{
    for(int i = 0; i < snl->n ; i++)
        var->Je[i] = evaluator_evaluate(snl->Bf[i], snl->n, snl->names, var->x0);
}

void substituteX(SistNl_t *snl, SnlVar_t *nt)
{
    double aux;
    char *func;
    // H[X]
    for(int i = 0; i < snl->n ; i++)
    {
        for(int j = 0; j < snl->n; j++){
            aux = evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, nt->x0);

            if(isnan(aux)) {
                // func = evaluator_get_string(snl->Hf[i][j]); // cursed
                // printf("func[%i][%i]: %s\n", i, j, func);
            }
            // printf("TA AKI O ERROOOOOOOOOO\n");
            nt->He[i][j] = aux;
        }

        // J[X]
        nt->Je[i] = evaluator_evaluate(snl->Bf[i], snl->n, snl->names, nt->x0);
    }

}

void snl2sl(SistNl_t *snl, SnlVar_t *nt)
{
    for(int i = 0; i < snl->n; ++i) {
        for(int j = 0; j < snl->n; ++j)
        {
            nt->sl->A[i][j] = nt->He[i][j];
        }
        nt->sl->b[i] = - nt->Je[i];
        // sl->t[i] ; // trocas salvas no pivo, dentro de nt->t[] (nao altera)
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

    /*
    // JACOBIANA
    printf("\nJACOBIANA\n#");
    for(int i = 0; i < S->n; i++)
        printf(" %f ", S->Je[i]);
    printf("\n");
    // for(int i = 0 ; i < snl->n; i++)
    //     printf("jacob %i : %f\n",i , evaluator_evaluate(snl->Bf[i], snl->n, names, values));

        
    // HESSIANA
    printf("\nHESSIANA\n#");
    for(int i = 0; i < S->n; i++){
        for(int j = 0; j < S->n; j++){
            printf("  %f  ", S->He[i][j]);
        }
        printf("\n");
    }
    // for(int i = 0 ; i < snl->n; i++){
    //     for(int j = 0; j < snl->n; j++)
    //         printf("%3f  ", evaluator_evaluate(snl->Hf[i][j], snl->n, snl->names, old_values));
    //     printf("\n");}
    */

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
    // func = evaluator_get_string(snl.Hf[0][0]); // cursed
    // printf("func[%i][%i]: %s\n", 0, 0, func);
    // for(int i = 0; i < snl.n; i++) printf("%s ", snl.names[i]); printf("\n");
        

    printf("-------------------------------------------------------------------------------\n");
}

void printCol(double* pto, int i, int max)
{
    if(i < max){
        if (isnan(pto[i]) || isinf(pto[i]))
            printf("%1.14e\t\t\t| ", pto[i]);
        else
            printf("%1.14e\t| ", pto[i]);
    }
    else
        printf("\t\t\t| ");

}
