touch newtonPC.c
make
make clear
echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

likwid-perfctr -O -C 3 -g L3 -m ./newtonPC < data/funcoesrosenbrock.dat > data/L3.log
likwid-perfctr -C 3 -g L3 -m ./newtonPC < data/funcoesrosenbrock.dat > data/L3.dat

echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
