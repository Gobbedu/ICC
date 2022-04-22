#!/bin/bash
touch newtonPC.c
make 
make clear

METRICA="L3 L2CACHE FLOPS_DP"
FUNCAO="../data/funcoesrosenbrock.dat"
TESTE="teste.dat"

RODAR=${FUNCAO}

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
# # TEMPO DE EXECUCAO
echo "calculando tempo padrao"
./newtonPC p ../data/csvs/OPT_tempoPadrao.csv  < ${RODAR}

echo "calculando tempo inexato"
./newtonPC i ../data/csvs/OPT_tempoInexato.csv < ${RODAR}



# # METRICAS LIKWID
for m in ${METRICA}
do
    echo "calculando likwid LOG ${m}"
    likwid-perfctr -O -C 3 -g ${m} -m ./newtonPC < ${RODAR} > ../data/logs/OPT_${m}.log
    
    # echo "calculando likwid DAT ${m}"
    # likwid-perfctr    -C 3 -g ${m} -m ./newtonPC < ${RODAR} > ../data/dats/OPT_${m}.dat
done
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

