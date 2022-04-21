/********************************************************
*   Resolução de Sistemas Nao Lineares
*   Eduardo Gobbo Willi GRR20203892 V.G. & Dante Eleuterio dos Santos GRR20206686
*   CI1164 - Introducao a Computacao Cientifica
*
*   ./newtonPC < funcoes.dat > saida_nossa.dat ou
*   ./newtonPC {_method(char) (i=inexato/ p=padrao)} {saida}
********************************************************/
#include "utils.h"
#include "NewtonPadrao.h"
#include "NewtonInexato.h"
// #define LIKWID_PERFMONI  -> Inclui likwid na compilacao, mudar em utils.h
// #define FULLPRINT_ON   //-> Imprime resultado de cara iteração
// #define _method 'o'// QUAL METODO EXECUTAR: p=newtonPadrao \ i=newtonInexato


int main(int argc, char **argv) {
    #ifdef LIKWID_PERFMONI
    LIKWID_MARKER_INIT;
    #endif

    char _method;

    SistNl_t *snl;
    Tempo_t tPadrao, tInexat;              // tempo de cada metodo

    initTempo(&tPadrao);
    initTempo(&tInexat);

    // tPadrao.derivadas = tInexat.derivadas = 0;
    // tPadrao.Gradiente = tInexat.Gradiente = 0;
    // tPadrao.Hessiana = tInexat.Hessiana = 0;
    // tPadrao.totalMetodo = tInexat.totalMetodo = 0;
    // tPadrao.totalSL = tInexat.totalSL = 0;

    double *respPadrao, *respInexat; 

    int iterPadrao = 0,    // por referencia, numerero de iteracoes (para imprimir)
        iterInexat = 0;

    FILE *saida;
    if(argc == 3)
    {
        saida = fopen(argv[2], "w");
        _method = argv[1][0];
    }    
    else    
        saida = stdout;

    // saida csv para 1 metodo separador ;
    fprintf(saida, "Aplicacao_metodo_Newton; Calculo_Gradiente; Calculo_Hessiana; Resolucao_Sistema_Linear\n");
    while(snl = lerSistNL())
    {   

        #ifdef FULLPRINT_ON
            fprintf(saida, "%i\n", snl->n);         // grau da funcao
            fprintf(saida, "%s\n", snl->funcao);    // a funcao
        #endif

        // calcula o He & o Ge dentro de cada metodo usando np/nm/ni
        if(_method == 'p'){
            respPadrao = malloc(sizeof(double) * snl->iteracao);
            NewtonPadrao(snl, respPadrao, &tPadrao, &iterPadrao);
        }
        else if(_method == 'i'){
            respInexat = malloc(sizeof(double) * snl->iteracao);
            NewtonInexato(snl,respInexat,&tInexat,&iterInexat);
        }
        else{
            respPadrao = malloc(sizeof(double) * snl->iteracao);
            NewtonPadrao(snl, respPadrao, &tPadrao, &iterPadrao);
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
            fprintf(saida, "%f; %f; %f; %f\n", 
            tPadrao.totalMetodo, tPadrao.Gradiente, tPadrao.Hessiana, tPadrao.totalSL);
            free(respPadrao);
        }
        else if(_method == 'i'){
            fprintf(saida, "%f; %f; %f; %f\n", 
            tInexat.totalMetodo, tInexat.Gradiente, tInexat.Hessiana, tInexat.totalSL);
            free(respInexat);
        }
        else{
            free(respPadrao);
            free(respInexat);
        }

        liberaSistNl(snl);
    }    
    if(argc == 3)
	    fclose(saida);

    #ifdef LIKWID_PERFMONI
    LIKWID_MARKER_CLOSE;
    #endif

    return 0;
}
