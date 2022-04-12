/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
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

    FILE *saida;
    if(argc == 3)   saida = fopen(argv[2], "w");

    while(snl = lerSistNL())
    {      
        snl->f = evaluator_create(snl->funcao);
        assert(snl->f);
        genNames(snl);          // nomes das variaveis x1, x2 .. xn

        // tPadrao.derivadas = timestamp();
        genGradiente(snl);   // derivadas de f c/ respeito a x1, x2 .. xn
        genHessiana(snl);       // possiveis combinacoes de segunda derivada
        // tPadrao.derivadas = timestamp() - tPadrao.derivadas;

        if(argc == 3){
            fprintf(saida, "%i\n", snl->n);
            fprintf(saida, "%s\n", snl->funcao);
        }
        else{
            printf("%i\n", snl->n);
            printf("%s\n", snl->funcao);
	    }

        respPadrao = malloc(sizeof(double) * snl->iteracao);
        respModifi = malloc(sizeof(double) * snl->iteracao);
        respInexat = malloc(sizeof(double) * snl->iteracao);

        // calcula o He & o Ge dentro de cada metodo usando np/nm/ni
        NewtonPadrao(snl, respPadrao, &tPadrao, &iterPadrao);
        NewtonModificado(snl, respModifi, &tModifi, &iterModifi);
        NewtonInexato(snl,respInexat,&tInexat,&iterInexat);

    	if(argc == 3)
            fprintf(saida, "#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");
	    else
	        printf("#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");

        // --------LOOP PRINT-------- //
        for(int i = 0; i < snl->iteracao; i++)
        {
            // imprime em -o <saida>
            if(argc == 3)
                fprintf(saida, "%-12d \t| ", i); // imprime iteração
            else
                printf("%-12d \t| ", i); // imprime iteração

            printCol(respPadrao, i, iterPadrao, argc, saida);
            printCol(respModifi, i, iterModifi, argc, saida);
            printCol(respInexat, i, iterInexat, argc, saida);

            if(argc == 3)
                fprintf(saida, "\n");
            else
                printf("\n");

            // se todos acabaram
            if( (i+1 >= iterPadrao) && (i+1 >= iterModifi) && (i+1>=iterInexat) )
                break;
        }

        if( argc == 3){
            fprintf(saida, "Tempo total \t| %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.totalMetodo, tModifi.totalMetodo, tInexat.totalMetodo);
            fprintf(saida, "Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.derivadas,tModifi.derivadas,tInexat.derivadas);
            fprintf(saida, "Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e  |\n#\n\n",tPadrao.totalSL,tModifi.totalSL,tInexat.totalSL);
        }
        else{
            printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.totalMetodo, tModifi.totalMetodo, tInexat.totalMetodo);
            printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.derivadas,tModifi.derivadas,tInexat.derivadas);
            printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e  |\n#\n\n",tPadrao.totalSL,tModifi.totalSL,tInexat.totalSL);
        }
        // LIBERA respMETODO
        free(respPadrao);
        free(respModifi);
        free(respInexat);

        // liberar matheval antes de destruir sistema
        liberaMatheval(snl);
        liberaSistNl(snl);
    }    
    
    if(argc == 3)   
        fclose(saida);

    return 0;
}
