import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

"""
    script para gerar plots em python
    tem como extrair informacao para plot de .csv com pandas
    
    especificacao .csv == separador ( ; ) ponto e virgula
    para cada csv gerado, utilizar estas colunas:
    Aplicacao_metodo_Newton; Calculo_Gradiente; Calculo_Hessiana; Resolucao_Sistema_Linear
    .                        .                  .                 .
    .                        .                  .                 .
    .                        .                  .                 .

    dataframe.plot() plota todas as colunas do csv
    dataframe["Nome Coluna"].plot() plota a Coluna do csv
    
    usar plt.scatter pra funcao que plota com pontos
    
    plt.show() mostra o grafico mas nao salva
    plt.savefig("nome_do_plot.png") salva o grafico em file
    
    rodar: $> python3 plotter.py

links pra ver como faz as parada:
pyplot[0]: https://matplotlib.org/3.1.1/tutorials/introductory/pyplot.html#sphx-glr-tutorials-introductory-pyplot-py
pandas[1]: https://pandas.pydata.org/docs/getting_started/intro_tutorials/04_plotting.html
"""

# diferentes dataframes
tempo_execucao        = pd.read_csv('data/auxiliar.csv')
banda_memoria         = pd.read_csv('data/auxiliar.csv')
cache_miss            = pd.read_csv('data/auxiliar.csv')
operacoes_aritmeticas = pd.read_csv('data/auxiliar.csv')



# fazer um grafico por vez
tempo_execucao["Pulse"].plot()                      # plota a coluna especificada
plt.yscale('log')                                   # especificado no trabalho
plt.title('teste com pyplot')                       # nome do grafico gerado, aparece na img dpois

plt.show()