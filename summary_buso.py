#!/usr/bin/env python3
# -*- coding: utf-8 -*- 
import glob,argparse,os
import pandas as pd
import numpy as np
#################################################
#北京组学生物科技有限公司
#python 语言入门与基础绘图：
# https://zzw.xet.tech/s/2bdd89
#生信基础课:
#https://zzw.xet.tech/s/1bWTkv
####################################################

parser = argparse.ArgumentParser(description='This script was used to stat gene busco summary result dir. ')
parser.add_argument('-d', '--dir', dest='dir', required=True, help='input  busco summary plot_dir dir')
parser.add_argument( '--run_type', dest='run_type', required=True, help='set busco run type, for example : specific')
parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default %(default)s')
parser.add_argument('-p','--prefix',dest='prefix',required=False,default="busco_summary",help='specify the output file prefix,default %(default)s')

args = parser.parse_args()
###################################################################

if not os.path.exists(args.outDir): os.mkdir(args.outDir)
args.outDir=os.path.abspath(args.outDir)



def load_data(plot_dir,run_type):
    """

    :return:
    """
    data={"Category":["Complete BUSCOs (C)","Complete and single-copy BUSCOs (S)","Complete and duplicated BUSCOs (D)","Fragmented BUSCOs (F)","Missing BUSCOs (M)","Total BUSCO groups searched"]}
    df=pd.DataFrame(data,index=np.arange(6))
#    print(df)
    for f in glob.glob("%s/short_summary.%s.*.*.txt" % (plot_dir, run_type)):
        try:
            
            content = open(f)
            comp = 0
            COMP= 0
            dupl = 0
            frag = 0
            miss = 0
            for line in content:
                if "Complete and single-copy BUSCOs" in line:
                    comp = int(line.split("\t")[1])
                elif "Complete and duplicated BUSCOs" in line:
                    dupl = int(line.split("\t")[1])
                elif "Complete BUSCOs" in line:
                    COMP = int(line.split("\t")[1])
                elif "Fragmented BUSCOs" in line:
                    frag = int(line.split("\t")[1])
                elif "Missing BUSCOs" in line:
                    miss = int(line.split("\t")[1])
                    
            species=f.split("/")[-1].split(".")[2]
            
            total = comp + dupl + frag + miss
            
            COMP_pc ="%s(%s%%)"%("{:,d}".format(COMP),"{:,.1f}".format(round(COMP / float(total) * 100, 1))) 
            comp_pc = "%s(%s%%)"%("{:,d}".format(comp),"{:,.1f}".format(round(comp / float(total) * 100, 1))) 
            dupl_pc = "%s(%s%%)"%("{:,d}".format(dupl),"{:,.1f}".format(round(dupl / float(total) * 100, 1))) 
            frag_pc = "%s(%s%%)"%("{:,d}".format(frag),"{:,.1f}".format(round(frag / float(total) * 100, 1))) 
            miss_pc = "%s(%s%%)"%("{:,d}".format(miss),"{:,.1f}".format(abs(round(100 - round(comp / float(total) * 100, 1) - round(dupl / float(total) * 100, 1) - round(frag / float(total) * 100, 1), 1)))) 
            
            print("Loaded %s successfully" % f)
        except IOError:
            print("Impossible to use the file %s" % f)
        
        df[species]=[COMP_pc,comp_pc, dupl_pc, frag_pc, miss_pc,"{:,d}".format(total)]
        
    return df

d=load_data(args.dir,args.run_type)
d.to_csv(args.outDir+"/"+args.prefix+".tsv",sep="\t",index=False)
