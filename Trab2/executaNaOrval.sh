touch newtonPC.c
make
make clear
echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
likwid-perfctr -C 3 -g L3 -m ./newtonPC < data/teste.dat
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
