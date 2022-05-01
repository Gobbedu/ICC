# Trabalho 1 ICC

### Feito por 
Eduardo Gobbo Willi Vasconcellos Gonçalves GRR20203892
&&
Dante Eleutério dos Santos GRR20206686


### Fluxo de Execução

Esse código busca resolver sistemas não lineares por meio de 3 variações do método iterativo de Newton: A primeira utiliza Eliminação de Gauss
para encontrar o x aproximado de cada iteração, a segunda utiliza o método da fatoração LU para tal e atualiza a matriz Hessiana a cada n iterações 
sendo n o número de variaveis e o terceiro método utiliza o método de Gauss Seidel para aproximar o x.
Os sistemas devem ser passados pela entrada padrão (stdin) ao chamar o executavel por meio de arquivo .dat e opcionalmente pode-se definir um arquivo para guardar a saída do código utilizando "-o" ou o ">" na chamada do executavel.
Ao fim do programa será imprimido os valores aproximados de f(x) a cada iteração de cada método e os tempos levados para rodar o método por completo, para calcular suas derivadas e para calcular o Sistema Linear.

Para reproduzir os resultados aqui obtidos, pode-se rodar o script full_script.sh.
Ajustes finos nos arquivos full_script.sh, extractinfo.py, plotter.py, utils.h e newtonPC.c podem ser necessários. 
Tais ajustes consistem em:
    - Especificar o diretório onde serão salvos o csv do likwid e a métrica de tempo em full_script.sh, como SAIDACSV e SAIDALOG
    - Especificar o diretórios de SAIDALOG em extractinfo.py na variavel logp, filtrando os dados do likwid em um csv utilizavel no diretório especificado por curp
    - Especificar em plotter.py:
        - o diretório onde se encontram os diretórios filtrados na variável src (metodo; grad; hess; sistlin)
        - o diretorio onde salvar os csvs que comparam o desempenho (otimizado; nao otimizado)
        - o local onde se salvar as imagens dos graficos na variavel out
    - Para que o script rode É NECESSÁRIO que a biblioteca likwid esteja instalada, e que seja possivel mudar o modo de operação da cpu com echo "perfomance/powersave" > scaling_governor. 

