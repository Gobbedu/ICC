#!/bin/bash
touch newtonPC.c
make
make clear

FUNCAO="data/funcoesrosenbrock.dat"
PREP="data/teste.dat"
RODAR=${PREP}

METRICA="L3 L2CACHE FLOPS_DP FLOPS_AVX"

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
for m in ${METRICA}
do
    likwid-perfctr -O -C 3 -g ${m} -m ./newtonPC < ${RODAR} > data/teste${m}.log
    likwid-perfctr    -C 3 -g ${m} -m ./newtonPC < ${RODAR} > data/teste${m}.dat
done
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

