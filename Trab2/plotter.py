from email.utils import collapse_rfc2231_value
from cv2 import rotate
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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

links pra ver como faz as parada:
pandas[1]: https://pandas.pydata.org/docs/getting_started/intro_tutorials/04_plotting.html

no final acho q nn vai precisar do pyplot
pyplot[0]: https://matplotlib.org/3.1.1/tutorials/introductory/pyplot.html#sphx-glr-tutorials-introductory-pyplot-py
"""


csv_file = 'data/timestampInexat.csv'
metodo = 'Inexato' # inexato ou padrao
saida_grafico = 'data/timestampInexat.png'
quero_ver = False    # se true mostra, se falso baixa em saida grafico

col1 = "Aplicacao_metodo_Newton"
col2 = "Calculo_Gradiente"
col3 = "Calculo_Hessiana"
col4 = "Resolucao_Sistema_Linear"

# diferentes dataframes
tempo_execucao = pd.read_csv(sep=';', filepath_or_buffer=csv_file)
# banda_memoria         = pd.read_csv('data/auxiliar.csv')
# cache_miss            = pd.read_csv('data/auxiliar.csv')
# operacoes_aritmeticas = pd.read_csv('data/auxiliar.csv')

# No eixo das abcissas os gráficos representam a dimensão N da Função de Rosenbrock
x = [10, 32, 50, 64, 100, 128, 200, 250, 256, 300, 400, 512, 600, 1000, 1024, 2000, 2048, 3000, 4096]
# fazer um grafico por vez

tempo_execucao.plot()

# o que vai ser medido
plt.ylabel("Tempo de execução")
# titulo do grafico
plt.title(f'Tempo Execucao Newton {metodo}')                       # nome do grafico gerado, aparece na img dpois

# NAO MUDA DAKI PRA BAIXO
plt.subplots_adjust(bottom=0.14)
plt.xticks(tempo_execucao.index, x[:len(tempo_execucao.index)], rotation=50) # magica em python  
plt.legend(["Total Metodo", "Gradiente", "Hessiana", "Sist. Linear"])
plt.xlabel("dimensão N da Função de Rosenbrock")
plt.yscale('log')                                               # especificado no Inexatotrabalho

# mostra grafico
if(quero_ver):
    plt.show()
else:
    # salvar grafico em:
    plt.savefig(saida_grafico)