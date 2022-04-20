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
#include <likwid.h>
// #define ROSENBROCK -> mudar em utils.h
// #define FULLPRINT_ON


int main(int argc, char **argv) {
    SistNl_t *snl;
    Tempo_t tPadrao, tInexat;              // tempo de cada metodo

    // initTempo(&tPadrao);
    // initTempo(&tInexat);

    // tPadrao.derivadas = tInexat.derivadas = 0;
    // tPadrao.Gradiente = tInexat.Gradiente = 0;
    // tPadrao.Hessiana = tInexat.Hessiana = 0;
    // tPadrao.totalMetodo = tInexat.totalMetodo = 0;
    // tPadrao.totalSL = tInexat.totalSL = 0;

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
    LIKWID_MARKER_INIT;
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

        char _method = 'p';

        // calcula o He & o Ge dentro de cada metodo usando np/nm/ni
        if(_method == 'p'){
            respPadrao = malloc(sizeof(double) * snl->iteracao);
            NewtonPadrao(snl, respPadrao, &tPadrao, &iterPadrao);
        }
        else if(_method == 'i'){
            respInexat = malloc(sizeof(double) * snl->iteracao);
            NewtonInexato(snl,respInexat,&tInexat,&iterInexat);
        }
        
        #ifdef FULLPRINT_ON
        fprintf(saida, "#Iteração \t| Newton Padrão \t| Newton Inexato\n");
        // --------LOOP PRINT-------- //
        for(int i = 0; i < snl->iteracao; i++)
        {
            fprintf(saida, "%-12d \t| ", i); // imprime iteração
            printCol(respPadrao, i, iterPadrao, saida);
            printCol(respInexat, i, iterInexat, saida);
  	        fprintf(saida, "\n");
            // se todos acabaram
            if( (i+1 >= iterPadrao) && (i+1 >= iterModifi) && (i+1>=iterInexat) )
                break;
        }
        // saida com os 2 metodos
        fprintf(saida, "Tempo total \t| %1.14e\t| %1.14e\t| \n",tPadrao.totalMetodo,  tInexat.totalMetodo);
        fprintf(saida, "Tempo derivadas | %1.14e\t| %1.14e\t| \n",tPadrao.derivadas, tInexat.derivadas);
        fprintf(saida, "Tempo SL \t| %1.14e\t| %1.14e\t| \n#\n\n",tPadrao.totalSL, tInexat.totalSL);
        #endif

        // saida para csv com metodo Padrao
        // LIBERA respMETODO
        if(_method =='p'){
            // fprintf(saida, "%f; %f; %f; %f\n", 
            // tPadrao.totalMetodo, tPadrao.Gradiente, tPadrao.Hessiana, tPadrao.totalSL);
            free(respPadrao);
        }
        else if(_method == 'i'){
            fprintf(saida, "%f; %f; %f; %f\n", 
            tInexat.totalMetodo, tInexat.Gradiente, tInexat.Hessiana, tInexat.totalSL);
            free(respInexat);
        }

        // liberar matheval antes de destruir sistema
        #ifndef ROSENBROCK
            liberaMatheval(snl);
        #endif

        liberaSistNl(snl);
    }    
    LIKWID_MARKER_CLOSE;
    if(argc == 3)
	    fclose(saida);
    return 0;
}
