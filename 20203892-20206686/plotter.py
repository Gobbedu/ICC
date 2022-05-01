#!/bin/python3
from os import minor
import matplotlib.pyplot as plt
import pandas as pd
import csv

"""
    script para gerar plots em python
    tem como extrair informacao para plot de .csv com pandas
    
    especificacao .csv == separador ( ; ) ponto e virgula
    para cada csv gerado, utilizar estas colunas:
    Aplicacao_metodo_Newton; Calculo_Gradiente; Calculo_Hessiana; Resolucao_Sistema_Linear
    TUDO                     .                  .                 .
    .                        .                  .                 .
    .                        .                  .                 .

    dataframe.plot() plota todas as colunas do csv
    dataframe["Nome Coluna"].plot() plota a Coluna do csv
    
    usar plt.scatter pra funcao que plota com pontos
    
    plt.show() mostra o grafico mas nao salva
    plt.savefig("nome_do_plot.png") salva o grafico em file
    
    rodar: $> python3 plotter.py
    rodar(todo): $> python3 plotter.py nomeMetodo nomeMedicao entrada.csv saida.csv

links pra ver como faz as parada:
pandas[1]: https://pandas.pydata.org/docs/getting_started/intro_tutorials/04_plotting.html
pyplot[0]: https://matplotlib.org/3.1.1/tutorials/introductory/pyplot.html#sphx-glr-tutorials-introductory-pyplot-py
"""

def plotter(input_file, out_plot, nameMetodo, metrica, log, save_plot):
    # saida_grafico = open(out_plot, 'w')
    input_csv = open(input_file, 'r')

    legenda=["Total Metodo", "Gradiente", "Hessiana", "Sist. Linear"]

    # diferentes dataframes
    df = pd.read_csv(sep=';', filepath_or_buffer=input_csv)

    # No eixo das abcissas os gráficos representam a dimensão N da Função de Rosenbrock
    x = [10, 32, 50, 64, 100, 128, 200, 250, 256, 300, 400, 512, 600, 1000, 1024, 2000, 2048, 3000, 4096]
    # fazer um grafico por vez

    # cores = list(islice(cycle(['b', 'g', 'orange', 'r']), None, len(df)))
    # df.plot(style='.-', color=cores)
    df.plot(style='.-')
    
    # o que vai ser medido
    plt.ylabel(metrica)
    # titulo do grafico
    plt.title(nameMetodo)    # nome do grafico gerado, aparece na img dpois

    # NAO MUDA DAKI PRA BAIXO
    plt.subplots_adjust(bottom=0.15, left=0.14, right=0.96)
    plt.xticks(df.index, x[:len(df.index)], rotation=50) # magica em python  
    # plt.legend(legenda)
    plt.grid()
    plt.xlabel("dimensão N da Função de Rosenbrock")
    if log:
        plt.yscale('log')                                               # especificado no Inexatotrabalho

    if (save_plot):
        plt.savefig(out_plot)  # salvar grafico em:
    else:
        plt.show()  # mostra grafico
        
    plt.close()
        

def showoffy(metrica):
    swicther = {
        'tempo': 'tempo em ms',
        'L3': 'L3 Memory bandwidth [MBytes/s]',
        'L2': 'L2 data cache miss ratio',
        'FLOPS_DP': 'FLOPS DP [MFLOP/s]',
        'FLOPS_AVX': 'FLOPS AVX [MFLOP/s]'
    }
    
    return swicther.get(metrica)


def showoffx(marker):
    swicther = {
        'METODO': 'Metodo',
        'GRAD': 'Gradiente',
        'HESS': 'Hessiana',
        'SISTLIN': 'Sistema Linear',
    }
    
    return 'Avalicação '+swicther.get(marker)+': '


    
def MABIKI(src, dest, metr_NoOpt_I, metr_Opt_I, metr_NoOpt_P, metr_Opt_P, MP, MI, GP, GI, HP, HI, SLP, SLI):
    aux1 = open(src+metr_NoOpt_I )
    NoOptI = csv.reader(aux1, delimiter=';')
    l_NoOptI = list(NoOptI)
    
    aux2 = open(src+metr_NoOpt_P )
    NoOptP = csv.reader(aux2, delimiter=';')
    l_NoOptP = list(NoOptP)
    
    aux3 = open(src+metr_Opt_I )
    OptI   = csv.reader(aux3, delimiter=';')
    l_OptI = list(OptI)
    
    aux4 = open(src+metr_Opt_P )
    OptP   = csv.reader(aux4, delimiter=';')
    l_OptP = list(OptP)

    # d = 'data/csvs/curated_csv/'
    fim = '.csv'
    
    auxMP_out = open(dest+MP+fim, "w")
    auxGP_out = open(dest+GP+fim, "w")
    auxHP_out = open(dest+HP+fim, "w")
    auxSLP_out= open(dest+SLP+fim, "w")
    
    auxMI_out = open(dest+MI+fim, "w")
    auxGI_out = open(dest+GI+fim, "w")
    auxHI_out = open(dest+HI+fim, "w")
    auxSLI_out= open(dest+SLI+fim, "w") 
    
    MP_out = csv.writer(auxMP_out,  delimiter=';')
    GP_out = csv.writer(auxGP_out,  delimiter=';')
    HP_out = csv.writer(auxHP_out,  delimiter=';')
    SLP_out = csv.writer(auxSLP_out, delimiter=';')

    MI_out = csv.writer(auxMI_out,  delimiter=';')
    GI_out = csv.writer(auxGI_out,  delimiter=';')
    HI_out = csv.writer(auxHI_out,  delimiter=';')
    SLI_out = csv.writer(auxSLI_out, delimiter=';')

    headerr = "Otimizado; Não Otimizado\n"

    auxMP_out.write(headerr)
    auxGP_out.write(headerr)
    auxHP_out.write(headerr)
    auxSLP_out.write(headerr)

    auxMI_out.write(headerr)
    auxGI_out.write(headerr)
    auxHI_out.write(headerr)
    auxSLI_out.write(headerr)

    
    for i in range(1, 20):
        MP_out  .writerow([l_OptP[i][0], l_NoOptP[i][0]])
        GP_out  .writerow([l_OptP[i][1], l_NoOptP[i][1]])
        HP_out  .writerow([l_OptP[i][2], l_NoOptP[i][2]])
        SLP_out .writerow([l_OptP[i][3], l_NoOptP[i][3]])
        
        MI_out  .writerow([l_OptI[i][0], l_NoOptI[i][0]])
        GI_out  .writerow([l_OptI[i][1], l_NoOptI[i][1]])
        HI_out  .writerow([l_OptI[i][2], l_NoOptI[i][2]])
        SLI_out .writerow([l_OptI[i][3], l_NoOptI[i][3]])
  
  
   
def MABIKI2(metr_NoOpt, metr_Opt, M, G, H, SL):
    aux1 = open(metr_NoOpt)
    NoOptI = csv.reader(aux1, delimiter=';')
    l_NoOpt = list(NoOptI)
    
    aux2 = open(metr_Opt)
    OptI   = csv.reader(aux2, delimiter=';')
    l_Opt = list(OptI)
    

    # d = 'data/csvs/curated_csv/'
    fim = '.csv'
    
    auxM_out = open(M+fim, "w")
    auxG_out = open(G+fim, "w")
    auxH_out = open(H+fim, "w")
    auxSL_out= open(SL+fim, "w")
    
    M_out = csv.writer(auxM_out,  delimiter=';')
    G_out = csv.writer(auxG_out,  delimiter=';')
    H_out = csv.writer(auxH_out,  delimiter=';')
    SL_out = csv.writer(auxSL_out, delimiter=';')

    headerr = "Otimizado; Não Otimizado\n"

    auxM_out.write(headerr)
    auxG_out.write(headerr)
    auxH_out.write(headerr)
    auxSL_out.write(headerr)

    
    for i in range(1, 20):
        M_out  .writerow([l_Opt[i][0], l_NoOpt[i][0]])
        G_out  .writerow([l_Opt[i][1], l_NoOpt[i][1]])
        H_out  .writerow([l_Opt[i][2], l_NoOpt[i][2]])
        SL_out .writerow([l_Opt[i][3], l_NoOpt[i][3]])


      
def moero_shinso_yo(src, dest, t,n,s,metrica):
# MABIKI(src, dest, 
#       metr_NoOpt_I, metr_Opt_I, 
#       metr_NoOpt_P, metr_Opt_P, 
#       MP, MI,
#       GP, GI,
#       HP, HI,
#       SLP, SLI):
    for m in metrica:
        print(m)
        MABIKI(src, dest,
                t[0]+m+n[0]+'.csv', t[1]+m+n[0]+'.csv', 
                t[0]+m+n[1]+'.csv', t[1]+m+n[1]+'.csv', 
                m+s[0]+n[1], m+s[0]+n[0],
                m+s[1]+n[1], m+s[1]+n[0],
                m+s[2]+n[1], m+s[2]+n[0],
                m+s[3]+n[1], m+s[3]+n[0])

def curate_raw_csv(src, dest, t,n,s,metrica):
    for newton in n:
        for m in metrica:
            print(m)
            MABIKI2(
                    src+t[0]+m+newton+'.csv', src+t[1]+m+newton+'.csv',
                    dest+m+s[0]+newton,
                    dest+m+s[1]+newton,
                    dest+m+s[2]+newton,
                    dest+m+s[3]+newton)

    

if __name__ == "__main__":
    
    t = ['noOPT_', 'OPT_']
    n = ['Inexato', 'Padrao']
    s = ['METODO', 'GRAD', 'HESS', 'SISTLIN']
    metrica = ['tempo', 'L3', 'L2', 'FLOPS_DP', 'FLOPS_AVX']
    
    src  = 'data/csvs/raw_csv/'
    curp = 'data/csvs/curated_csv/'
    # out  = 'data/plots/opt_plots/'
    out = 'data/plots/'
    
    # parse raw csv to compare opt to noopt
    # moero_shinso_yo(src, curp, t,n,s,metrica)
    curate_raw_csv(src, curp, t,n,s,metrica)
    
    # estrutura files [metrica][marker][metodo]
    for m in metrica:
        for marker in s:
            for metodo in n:
                plotter(curp+m+marker+metodo+'.csv',
                        out+m+'_'+marker+metodo+'.png',
                        showoffx(marker)+' Newton ' +metodo,
                        showoffy(m), 
                        m == 'tempo',
                        True
                        )
    

  
