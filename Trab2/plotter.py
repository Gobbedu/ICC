#!/bin/python3
from itertools import cycle, islice
import matplotlib.pyplot as plt
import pandas as pd
from sqlalchemy import false

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

def plotter(input_csv, legenda, out_plot, nameMetodo, metrica, log, save_plot):
    saida_grafico = out_plot

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
    plt.title(f'{nameMetodo}')    # nome do grafico gerado, aparece na img dpois

    # NAO MUDA DAKI PRA BAIXO
    plt.subplots_adjust(bottom=0.14)
    plt.xticks(df.index, x[:len(df.index)], rotation=50) # magica em python  
    plt.legend(legenda)
    plt.xlabel("dimensão N da Função de Rosenbrock")
    if log:
        plt.yscale('log')                                               # especificado no Inexatotrabalho

    if (save_plot):
        plt.savefig(saida_grafico)  # salvar grafico em:
    else:
        plt.show()  # mostra grafico
        
        
    
if __name__ == "__main__":
    # metrica eh tempo
    TnoOptP = open('data/csvs/noOPT_tempoPadrao.csv', 'r')
    TnoOptI = open('data/csvs/noOPT_tempoInexato.csv', 'r')
    TOptP   = open('data/csvs/OPT_tempoPadrao.csv', 'r')
    TOptI   = open('data/csvs/OPT_tempoInexato.csv', 'r')
    
    # metrica eh L3
    L3noOptP = open('data/csvs/noOPT_L3Padrao.csv', 'r')
    L3noOptI = open('data/csvs/noOPT_L3Inexato.csv', 'r')
    L3OptP   = open('data/csvs/OPT_L3Inexato.csv', 'r')
    L3OptI   = open('data/csvs/OPT_L3Inexato.csv', 'r')
    
    # metrica eh L2
    L2P  = open('data/csvs/noOPT_L2Padrao.csv', 'r')
    L2I  = open('data/csvs/noOPT_L2Inexato.csv', 'r')
    L2oP = open('data/csvs/OPT_L2Inexato.csv', 'r')
    L2oI = open('data/csvs/OPT_L2Inexato.csv', 'r')
    
    
    
    legenda = ["Total Metodo", "Gradiente", "Hessiana", "Sist. Linear"]

    salve = True
    log = True
    metric = 'Tempo de Execução'
    plotter(TnoOptP, legenda, 'data/plots/noOpt_tempoPadrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(TnoOptI, legenda, 'data/plots/noOpt_tempoInexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(TOptP, legenda, 'data/plots/Opt_tempoPadrao.png.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(TOptI, legenda, 'data/plots/Opt_tempoInexato.png.png', 'Newton Inexato Otimizado', metric, log, salve)

    metric = 'Memory bandwidth [MBytes/s]'
    plotter(L3noOptP, legenda, 'data/plots/noOpt_L3Padrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(L3noOptI, legenda, 'data/plots/noOpt_L3Inexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(L3OptP, legenda, 'data/plots/Opt_L3Padrao.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(L3OptI, legenda, 'data/plots/Opt_L3Inexato.png', 'Newton Inexato Otimizado', metric, log, salve)

    log = False
    metric = 'data cache miss ratio'
    plotter(L2P, legenda, 'data/plots/noOpt_L2Padrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(L2I, legenda, 'data/plots/noOpt_L2Inexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(L2oP, legenda, 'data/plots/Opt_L2Padrao.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(L2oI, legenda, 'data/plots/Opt_L2Inexato.png', 'Newton Inexato Otimizado', metric, log, salve)

