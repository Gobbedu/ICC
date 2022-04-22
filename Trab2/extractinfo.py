#!/bin/python3

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

def parse_LOG(csv_file, padrao_file, inexat_file, line_of_interest):
    padrao_out = open(padrao_file, "w")
    inexat_out = open(inexat_file, "w")
    csv_in     = open(csv_file, "r")
    
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
    parse_LOG(csv_in, padrao_out, inexat_out, 10)
                    
def parse_L2(csv_in, padrao_out, inexat_out):
    parse_LOG(csv_in, padrao_out, inexat_out, 8)

def parse_FLOPS(csv_file, padraoDP_out, inexatDP_out, padraoAVX_out=None, inexatAVX_out=None):
    parse_LOG(csv_file, padraoDP_out, inexatDP_out, 6)
    parse_LOG(csv_file, padraoAVX_out, inexatAVX_out, 7)
    


if __name__ == "__main__":
    # data de nao otimizado
    raw_DP       = 'data/logs/noOPT_FLOPS_DP.log'
    raw_L3       = 'data/logs/noOPT_L3.log'
    raw_L2       = 'data/logs/noOPT_L2CACHE.log'
    L3padrao_csv = 'data/csvs/noOPT_L3Padrao.csv'
    L3inexato_csv= 'data/csvs/noOPT_L3Inexato.csv'
    L2padrao_csv = 'data/csvs/noOPT_L2Padrao.csv'
    L2inexato_csv= 'data/csvs/noOPT_L2Inexato.csv'
    DPP          = 'data/csvs/noOPT_FLOPS_DPPadrao.csv'
    DPI          = 'data/csvs/noOPT_FLOPS_DPInexato.csv'
    AVXP         = 'data/csvs/noOPT_FLOPS_AVX_Padrao.csv'
    AVXI         = 'data/csvs/noOPT_FLOPS_AVX_Inexato.csv'

    # data de otimizado
    raw_DPO = 'data/logs/OPT_FLOPS_DP.log'
    raw_L3O = 'data/logs/OPT_L3.log'
    raw_L2O = 'data/logs/OPT_L2CACHE.log'
    L3P     = 'data/csvs/OPT_L3Padrao.csv'
    L3I     = 'data/csvs/OPT_L3Inexato.csv'
    L2P     = 'data/csvs/OPT_L2Padrao.csv'
    L2I     = 'data/csvs/OPT_L2Inexato.csv'
    DPoP    = 'data/csvs/OPT_FLOPS_DPPadrao.csv'
    DPoI    = 'data/csvs/OPT_FLOPS_DPInexato.csv'
    AVXoP   = 'data/csvs/OPT_FLOPS_AVX_Padrao.csv'
    AVXoI   = 'data/csvs/OPT_FLOPS_AVX_Inexato.csv'
    
    # gera csv para plotter
    parse_L3(raw_L3, L3padrao_csv, L3inexato_csv)
    parse_L2(raw_L2, L2padrao_csv, L2inexato_csv)
    parse_FLOPS(raw_DP, DPP, DPI, AVXP, AVXI)

    parse_L3(raw_L3O, L3P, L3I)
    parse_L2(raw_L2O, L2P, L2I)
    parse_FLOPS(raw_DPO, DPoP, DPoI, AVXoP, AVXoI)

