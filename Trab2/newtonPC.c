/********************************************************
*    Resolução de Sistemas Nao Lineares
*    Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*    CI1164 - Introducao a Computacao Cientifica
*
*    ./newtonPC < funcoes.dat > saida_nossa.dat
********************************************************/
#include "utils.h"
// #include "NewtonModificado.h"
#include "NewtonPadrao.h"
#include "NewtonInexato.h"

// #define ROSENBROCK -> mudar em utils.h
// #define FULLPRINT_ON


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
    if(argc == 3)    
        saida = fopen(argv[2], "w");
    else    
        saida = stdout;

    // saida csv para 1 metodo separador ;
    fprintf(saida, "Aplicacao_metodo_Newton; Calculo_Gradiente; Calculo_Hessiana; Resolucao_Sistema_Linear\n");

    while(snl = lerSistNL())
    {   
        #ifndef ROSENBROCK   
            snl->f = evaluator_create(snl->funcao);
            assert(snl->f);
            genNames(snl);          // nomes das variaveis x1, x2 .. xn

            genGradiente(snl);   // derivadas de f c/ respeito a x1, x2 .. xn
            genHessiana(snl);       // possiveis combinacoes de segunda derivada
        #endif

        #ifdef FULLPRINT_ON
            fprintf(saida, "%i\n", snl->n);         // grau da funcao
            fprintf(saida, "%s\n", snl->funcao);    // a funcao
        #endif

        respPadrao = malloc(sizeof(double) * snl->iteracao);
        // respInexat = malloc(sizeof(double) * snl->iteracao);
        // respModifi = malloc(sizeof(double) * snl->iteracao);

        // calcula o He & o Ge dentro de cada metodo usando np/nm/ni
        // NewtonModificado(snl, respModifi, &tModifi, &iterModifi);
        NewtonPadrao(snl, respPadrao, &tPadrao, &iterPadrao);
        // NewtonInexato(snl,respInexat,&tInexat,&iterInexat);

        
        // fprintf(saida, "#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");
        // fprintf(saida, "#Iteração \t| Newton Padrão \t| Newton Inexato\n");


        #ifdef FULLPRINT_ON
        // --------LOOP PRINT-------- //
        for(int i = 0; i < snl->iteracao; i++)
        {
	    // imprime em -o <saida>
            fprintf(saida, "%-12d \t| ", i); // imprime iteração

            printCol(respPadrao, i, iterPadrao, saida);
            printCol(respInexat, i, iterInexat, saida);
            // printCol(respModifi, i, iterModifi, argc, saida);

  	        fprintf(saida, "\n");

            // se todos acabaram
            if( (i+1 >= iterPadrao) && (i+1 >= iterModifi) && (i+1>=iterInexat) )
                break;
        }
        #endif
        // saida com os 3 metodos
        // fprintf(saida, "Tempo total \t| %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.totalMetodo, tModifi.totalMetodo, tInexat.totalMetodo);
        // fprintf(saida, "Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e  |\n",tPadrao.derivadas,tModifi.derivadas,tInexat.derivadas);
        // fprintf(saida, "Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e  |\n#\n\n",tPadrao.totalSL,tModifi.totalSL,tInexat.totalSL);
        // saida com os 2 metodos
        // fprintf(saida, "Tempo total \t| %1.14e\t| %1.14e\t| \n",tPadrao.totalMetodo,  tInexat.totalMetodo);
        // fprintf(saida, "Tempo derivadas | %1.14e\t| %1.14e\t| \n",tPadrao.derivadas, tInexat.derivadas);
        // fprintf(saida, "Tempo SL \t| %1.14e\t| %1.14e\t| \n#\n\n",tPadrao.totalSL, tInexat.totalSL);

        // saida para csv com 1 metodo
        fprintf(saida, "%f; %f; %f; %f\n", 
        tPadrao.totalMetodo, tPadrao.Gradiente, tPadrao.Hessiana, tPadrao.totalSL);

        // LIBERA respMETODO
        free(respPadrao);
        // free(respModifi);
        // free(respInexat);

        // liberar matheval antes de destruir sistema
        #ifndef ROSENBROCK
            liberaMatheval(snl);
        #endif

        liberaSistNl(snl);
    }    

    if(argc == 3)
	    fclose(saida);
    return 0;
}
