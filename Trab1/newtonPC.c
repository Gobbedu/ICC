/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi V.G. & Dante Eleuterio dos Santos
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat
********************************************************/

#include "utils.h"
#include "NewtonModificado.h"
#include "NewtonPadrao.h"
#include "NewtonInexato.h"

// #define MAXIT_REFINAMENTO 5

// #define DEBUG_FLAG


int main(int argc, char **argv) {
    SistNl_t *snl;
    Tempo_t tPadrao, tModifi, tInexat;              // tempo de cada metodo

    initTempo(&tPadrao);
    initTempo(&tModifi);
    initTempo(&tInexat);

    double *respPadrao, *respModifi, *respInexat; 

    int iterPadrao = 0,    // por referencia, numerero de iteracoes (para imprimir)
        iterModifi = 0,
        iterInexat = 0;

    while(snl = lerSistNL())
    {      
        snl->f = evaluator_create(snl->funcao);
        assert(snl->f);
        genNames(snl);          // nomes das variaveis x1, x2 .. xn

        // tPadrao.derivadas = timestamp();
        genGradiente(snl);   // derivadas de f c/ respeito a x1, x2 .. xn
        genHessiana(snl);       // possiveis combinacoes de segunda derivada
        // tPadrao.derivadas = timestamp() - tPadrao.derivadas;

        printf("%i\n", snl->n);
        printf("%s\n", snl->funcao);

        respPadrao = malloc(sizeof(double) * snl->iteracao);
        respModifi = malloc(sizeof(double) * snl->iteracao);
        respInexat = malloc(sizeof(double) * snl->iteracao);
        // snlinfo(snl);

        // snl precisa de copia, muda o He & o Ge
        // calcula o He & o Ge dentro de cada metodo usando np/nm/ni

        NewtonPadrao(snl, respPadrao, &tPadrao, &iterPadrao);
        NewtonModificado(snl, respModifi, &tModifi, &iterModifi);
        NewtonInexato(snl,respInexat,&tInexat,&iterInexat);

        printf("#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");

        // --------LOOP PRINT-------- //
        for(int i = 0; i < snl->iteracao; i++)
        {
            printf("%-12d \t| ", i); // imprime iteração
            printCol(respPadrao, i, iterPadrao);
            printCol(respModifi, i, iterModifi);
            printCol(respInexat, i, iterInexat);
            printf("\n");

            // se todos acabaram
            if( (i+1 >= iterPadrao) && (i+1 >= iterModifi) && (i+1>=iterInexat) )
                break;
        }
        printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.totalMetodo, tModifi.totalMetodo, tInexat.totalMetodo);
        printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.derivadas,tModifi.derivadas,tInexat.derivadas);
        printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e  |\n#\n\n",tPadrao.totalSL,tModifi.totalSL,tInexat.totalSL);

        // LIBERA respMETODO
        free(respPadrao);
        free(respModifi);
        free(respInexat);

        // liberar matheval antes de destruir sistema
        liberaMatheval(snl);
        liberaSistNl(snl);
    }    
    return 0;
}
