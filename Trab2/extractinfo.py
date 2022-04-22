#!/bin/python3

import matplotlib.pyplot as plt
from numpy import pad
import pandas as pd

"""
    script para extrair metricas do .log
    
    formato:
    recebe .log de metrica
    
    devolve:
    sobre uma metrica (L3/L2/FLOPS)
               |  csv_P  |  csv_I  |
    -----------+---------+---------+
    Tmetodo    |    ~    |    ~    |
    Gradiente  |    ~    |    ~    |
    Hessiana   |    ~    |    ~    |
    Sist lin   |    ~    |    ~    |
    -----------+---------+---------+
    
""" 

def parse_LOG(csv_in, padrao_out, inexat_out, line_of_interest):
    padrao_out.write("METODO; GRAD; HESS; SISTLIN\n")
    inexat_out.write("METODO; GRAD; HESS; SISTLIN\n")
    lines = csv_in.readlines()

    metricsP = {"METODO":0,"GRAD":0, "HESS":0, "SISTLIN":0}
    metricsI = {"METODO":0,"GRAD":0, "HESS":0, "SISTLIN":0}
    keys = ["METODO","GRAD","HESS","SISTLIN"]
    for pos, line in enumerate(lines):
        if "Group 1 Metric" in line and "Raw" not in line:
            # for metric in metrics:
            # print(line, lines[pos+line_of_interest])
# """
            Nth = line.split(',')[1].split(' ')[1].split('_')
            value = lines[pos+line_of_interest].split(',')[1]
            # print(Nth, ":  ", value)
            for key in metricsI.keys():
                if key in Nth[0] and "Padrao" in Nth[0]:
                    metricsP[key] = value
                if key in Nth[0] and "Inexato" in Nth[0]:
                    metricsI[key] = value

            # se preencheu dict PADRAO, write && zera
            if all(val != 0 for val in metricsP.values()):   
                # print("METRICA p: ",metricsP.values())
                m = metricsP.get(keys[0])
                g = metricsP.get(keys[1])
                h = metricsP.get(keys[2])
                s = metricsP.get(keys[3])
                padrao_out.write(f"{m}; {g}; {h}; {s}\n".format(m, g, h, s))

                metricsP = dict.fromkeys(keys, 0)
                # print("padrao:", metricsP, end='\n')

            # se preencheu dict INEXATO, write && zera
            if all(val != 0 for val in metricsI.values()):
                # print("METRICA i: ",metricsP.values())
                m = metricsI.get(keys[0])
                g = metricsI.get(keys[1])
                h = metricsI.get(keys[2])
                s = metricsI.get(keys[3])
                inexat_out.write(f"{m}; {g}; {h}; {s}\n".format(m, g, h, s))
                metricsI = dict.fromkeys(keys, 0)
                # print("inexato:", metricsI, end='\n')
# """

def parse_L3(csv_in, padrao_out, inexat_out):
    parse_LOG(csv_in, padrao_out, inexat_out, 6)
                    
def parse_L2(csv_in, padrao_out, inexat_out):
    parse_LOG(csv_in, padrao_out, inexat_out, 8)

def parse_FLOPS(csv_file, padraoDP_out, inexatDP_out, padraoAVX_out=None, inexatAVX_out=None):
    # csv_in = open(csv_file, "r")
    parse_LOG(csv_file, padraoDP_out, inexatDP_out, 6)
    # csv_in = open(csv_file, "r")
    # parse_LOG(csv_in, padraoAVX_out, inexatAVX_out, 7)
    


if __name__ == "__main__":
    # data de nao otimizado
    raw_L3 = open("data/logs/noOPT_L3.log", "r")
    raw_L2 = open("data/logs/noOPT_L2CACHE.log", "r")
    raw_DP = open("data/logs/noOPT_FLOPS_DP.log")
    L3padrao_csv = open("data/csvs/noOPT_L3Padrao.csv", "w")
    L3inexato_csv= open("data/csvs/noOPT_L3Inexato.csv", "w")
    L2padrao_csv = open("data/csvs/noOPT_L2Padrao.csv", "w")
    L2inexato_csv= open("data/csvs/noOPT_L2Inexato.csv", "w")
    DPpadrao_csv = open("data/csvs/noOPT_FLOPS_DPPadrao.csv", "w")
    DPinexato_csv= open("data/csvs/noOPT_FLOPS_DPInexato.csv", "w")
    # AVXpadrao_csv= open("data/csvs/noOPTPadrao_FLOPS_AVX.csv", "w")
    # AVXinexat_csv= open("data/csvs/noOPTInexat_FLOPS_AVX.csv", "w")

    # data de otimizado
    raw_L3O = open('data/logs/OPT_L3.log')
    raw_L2O = open('data/logs/OPT_L2CACHE.log')
    raw_DPO = open('data/logs/OPT_FLOPS_DP.log')
    L3P = open("data/csvs/OPT_L3Padrao.csv", "w")
    L3I = open("data/csvs/OPT_L3Inexato.csv", "w")
    L2P = open("data/csvs/OPT_L2Padrao.csv", "w")
    L2I = open("data/csvs/OPT_L2Inexato.csv", "w")
    DPoP= open("data/csvs/OPT_FLOPS_DPPadrao.csv", "w")
    DPoI= open("data/csvs/OPT_FLOPS_DPInexato.csv", "w")
    # AVXpadrao_csv= open("data/csvs/noOPTPadrao_FLOPS_AVX.csv", "w")
    # AVXinexat_csv= open("data/csvs/noOPTInexat_FLOPS_AVX.csv", "w")
    
    # gera csv para plotter
    parse_L3(raw_L3, L3padrao_csv, L3inexato_csv)
    parse_L2(raw_L2, L2padrao_csv, L2inexato_csv)
    parse_FLOPS(raw_DP, DPpadrao_csv, DPinexato_csv)

    parse_L3(raw_L3O, L3P, L3I)
    parse_L2(raw_L2O, L2P, L2I)
    parse_FLOPS(raw_DPO, DPoP, DPoI)

