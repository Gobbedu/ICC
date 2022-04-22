#!/bin/python3
import matplotlib.pyplot as plt
import pandas as pd

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

def plotter(input_file, out_plot, nameMetodo, metrica, log, save_plot, legenda=None):
    # saida_grafico = open(out_plot, 'w')
    input_csv = open(input_file, 'r')

    if not legenda:
        legenda = ["Total Metodo", "Gradiente", "Hessiana", "Sist. Linear"]

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
        plt.savefig(out_plot)  # salvar grafico em:
    else:
        plt.show()  # mostra grafico
        
        
    
if __name__ == "__main__":
    # metrica eh tempo
    TnoOptP = 'data/csvs/noOPT_tempoPadrao.csv'
    TnoOptI = 'data/csvs/noOPT_tempoInexato.csv'
    TOptP   = 'data/csvs/OPT_tempoPadrao.csv'
    TOptI   = 'data/csvs/OPT_tempoInexato.csv'
    
    # metrica eh L3
    L3noOptP = 'data/csvs/noOPT_L3Padrao.csv'
    L3noOptI = 'data/csvs/noOPT_L3Inexato.csv'
    L3OptP   = 'data/csvs/OPT_L3Inexato.csv'
    L3OptI   = 'data/csvs/OPT_L3Inexato.csv'
    
    # metrica eh L2
    L2P  = 'data/csvs/noOPT_L2Padrao.csv'
    L2I  = 'data/csvs/noOPT_L2Inexato.csv'
    L2oP = 'data/csvs/OPT_L2Inexato.csv'
    L2oI = 'data/csvs/OPT_L2Inexato.csv'
    
    # metrica eh FLOPS DP
    FDPp  = 'data/csvs/noOPT_FLOPS_DPPadrao.csv'
    FDPi  = 'data/csvs/noOPT_FLOPS_DPInexato.csv'
    FDPop = 'data/csvs/OPT_FLOPS_DPPadrao.csv'
    FDPoi = 'data/csvs/OPT_FLOPS_DPInexato.csv'
    
    # FLOPS AVX
    FAVXp = 'data/csvs/noOPT_FLOPS_AVX_Padrao.csv'
    FAVXi = 'data/csvs/noOPT_FLOPS_AVX_Inexato.csv'
    FAVXop= 'data/csvs/OPT_FLOPS_AVX_Padrao.csv'
    FAVXoi= 'data/csvs/OPT_FLOPS_AVX_Inexato.csv'
    
    salve = True
    log = True
    metric = 'Tempo de Execução'
    plotter(TnoOptP, 'data/plots/noOpt_tempoPadrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(TnoOptI, 'data/plots/noOpt_tempoInexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(TOptP, 'data/plots/Opt_tempoPadrao.png.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(TOptI, 'data/plots/Opt_tempoInexato.png.png', 'Newton Inexato Otimizado', metric, log, salve)

    metric = 'Memory bandwidth [MBytes/s]'
    plotter(L3noOptP, 'data/plots/noOpt_L3Padrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(L3noOptI, 'data/plots/noOpt_L3Inexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(L3OptP, 'data/plots/Opt_L3Padrao.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(L3OptI, 'data/plots/Opt_L3Inexato.png', 'Newton Inexato Otimizado', metric, log, salve)

    log = False
    metric = 'data cache miss ratio'
    plotter(L2P, 'data/plots/noOpt_L2Padrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(L2I, 'data/plots/noOpt_L2Inexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(L2oP, 'data/plots/Opt_L2Padrao.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(L2oI, 'data/plots/Opt_L2Inexato.png', 'Newton Inexato Otimizado', metric, log, salve)

    log = True
    metric = 'FLOPS DP'
    plotter(FDPp, 'data/plots/noOpt_DP_Padrao.png', 'Newton Padrao não Otimizado', metric, log, salve)
    plotter(FDPi, 'data/plots/noOpt_DP_Inexato.png', 'Newton Inexato não Otimizado', metric, log, salve)
    plotter(FDPop, 'data/plots/Opt_DP_Padrao.png', 'Newton Padrao Otimizado', metric, log, salve)
    plotter(FDPoi, 'data/plots/Opt_DP_Inexato.png', 'Newton Inexato Otimizado', metric, log, salve)

    metric = 'FLOPS AVX'
    leg = ["Total Metodo", "Sist Linear"]
    plotter(FAVXp, 'data/plots/noOpt_AVX_Padrao.png', 'Newton Padrao não Otimizado', metric, log, salve, leg)
    plotter(FAVXi, 'data/plots/noOpt_AVX_Inexato.png', 'Newton Inexato não Otimizado', metric, log, salve, leg)
    plotter(FAVXop, 'data/plots/Opt_AVX_Padrao.png', 'Newton Padrao Otimizado', metric, log, salve, leg)
    plotter(FAVXoi, 'data/plots/Opt_AVX_Inexato.png', 'Newton Inexato Otimizado', metric, log, salve, leg)
