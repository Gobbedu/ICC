#!/bin/bash

touch newtonPC.C
make local
make clear

TESTE="data/teste.dat"
ROSEN="data/funcoesrosenbrock.dat"

ENTRADA=${ROSEN}

echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
./newtonPC p data/csvs/noOPT_tempoPadrao.csv  < ${ENTRADA}
./newtonPC i data/csvs/noOPT_tempoInexato.csv < ${ENTRADA}
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

