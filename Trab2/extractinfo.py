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
        if "Group 1 Metric" in line:
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


def parse_L3(log_file, padrao_out, inexat_out):
    parse_LOG(log_file, padrao_out, inexat_out, 10)
            
                    
def parse_L2(log_file, padrao_out, inexat_out):
    parse_LOG(log_file, padrao_out, inexat_out, 8)


def parse_FLOPS(log_file, padraoDP_out, inexatDP_out, padraoAVX_out=None, inexatAVX_out=None):
    parse_LOG(log_file, padraoDP_out, inexatDP_out, 6)
    parse_LOG(log_file, padraoAVX_out, inexatAVX_out, 7)
    


if __name__ == "__main__":
    logp = 'data/logs/opt_log/'
    curp = 'data/csvs/opt_csv/'
    
    tipos = ['noOPT_', 'OPT_']
    newt  = ['Padrao', 'Inexato']
    metrs = ['FLOPS_DP', 'L3', 'L2CACHE']
    
    # PADRAO => {TIPO}{METRICA}{METODO}{EXT}
    # ex: noOPT_L3.log -> noOPT_L3Padrao.csv,  noOPT_L3Inexato.csv

    for tipo in tipos:
        parse_L3(logp+tipo+'L3.log',
                 curp+tipo+'L3'+newt[0]+'.csv',
                 curp+tipo+'L3'+newt[1]+'.csv')
        
        parse_L2(logp+tipo+'L2CACHE.log',
                 curp+tipo+'L2'+newt[0]+'.csv',
                 curp+tipo+'L2'+newt[1]+'.csv')
        
        parse_FLOPS(logp+tipo+'FLOPS_DP.log', 
                    curp+tipo+'FLOPS_DP'+newt[0]+'.csv', 
                    curp+tipo+'FLOPS_DP'+newt[1]+'.csv',
                    curp+tipo+'FLOPS_AVX'+newt[0]+'.csv',
                    curp+tipo+'FLOPS_AVX'+newt[1]+'.csv')

