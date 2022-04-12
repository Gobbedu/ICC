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
