#!/bin/bash

METRICA="L3 L2CACHE FLOPS_DP"
FUNCAO="../data/funcoesrosenbrock.dat"
TESTE="../data/teste.dat"

SAIDACSV="../data/csvs/opt_csv/"
SAIDALOG="../data/logs/opt_log/"

# for python script
# raw log
# curated log




RODAR=${TESTE}

# ============== EXECUTA NAO OTIMIZADO ==============
cd Rosbrk/
pwd
touch newtonPC.c
make
make clear

# echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

# TEMPO DE EXECUCAO
echo "calculando tempo padrao NAO OTIMIZADO"
./newtonPC p ${SAIDACSV}noOPT_tempoPadrao.csv  < ${RODAR}

echo "calculando tempo inexato NAO OTIMIZADO"
./newtonPC i ${SAIDACSV}noOPT_tempoInexato.csv < ${RODAR}

# METRICAS LIKWID
for m in ${METRICA}
do
    echo "calculando likwid NAO OTIMIZADO ${m} ..."
    likwid-perfctr -O -C 3 -g ${m} -m ./newtonPC < ${RODAR} > ${SAIDALOG}noOPT_${m}.log
done

# echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

make purge
cd ..

# ============= EXECUTA OTIMIZADO =============
cd Rosbrk_Opt
pwd
touch newtonPC.C
make
make clear

# echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

# TEMPO DE EXECUCAO
echo "calculando tempo padrao OTIMIZADO "
./newtonPC p ${SAIDACSV}noOPT_tempoPadrao.csv  < ${RODAR}

echo "calculando tempo inexato OTIMIZADO "
./newtonPC i ${SAIDACSV}noOPT_tempoInexato.csv < ${RODAR}

# METRICAS LIKWID
for m in ${METRICA}
do
    echo "calculando likwid OTIMIZADO ${m} ..."
    likwid-perfctr -O -C 3 -g ${m} -m ./newtonPC < ${RODAR} > ${SAIDALOG}noOPT_${m}.log
done

# echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

make purge
cd ..

# =================== CRIA PLOTS ===================
pwd
./extractinfo.py
./plotter.py

