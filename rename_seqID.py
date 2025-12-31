# -*- coding: utf-8 -*- 
'''
This script was used to convert xml config file to html web page.

Author:
       gejy
Version:
        1.0;2022-7-25

python 语言入门与基础绘图：
 https://zzw.xet.tech/s/2bdd89
生信基础课:
https://zzw.xet.tech/s/1bWTkv

'''



############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2022.07.10
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.h5.xeknow.com/s/2G8tHr
###############################################################################################

import argparse
import os
import re
parser = argparse.ArgumentParser(description='This script is used to rename fasta and apg files')

parser.add_argument('-r','--rename_table',help='Please input rename file',required=True)
parser.add_argument('-f','--fasta',help='Please input fasta file',required=False)
parser.add_argument('-a','--agp',help='Please input file',required=False)
parser.add_argument('--remove_desc',help='remove description info from fasta id ',action='store_true',required=False)
#parser.add_argument('-t','--table',type=int,default=1,help=' genetic code :https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi, default table id 1',required=False)
parser.add_argument('-o','--out_dir',help='Please input complete out_put directory path',default = os.getcwd(),required=False)
parser.add_argument('-p','--prefix',default ='pep',required=False,help='Please specify the output, pep')
################################################################################
args = parser.parse_args()
dout=''
if os.path.exists(args.out_dir):
    dout=os.path.abspath(args.out_dir)
else:
    os.mkdir(args.out_dir)
    dout=os.path.abspath(args.out_dir)

#把chr_rename提取出来，第一列作为key，第二列作为value，得到一个字典name
f1 = open(args.rename_table, "r")
name={}
for line in f1:
    line=line.strip()
    if line == "":continue
    list1 = re.split('\s+',line)
    name[list1[0]] = list1[1]

if args.fasta:
#对fasta的每一行，去掉">"，比一下序列id在不在字典中
    f2 = open(args.fasta, "r")
    f4 = open(dout+"/"+args.prefix+".fa", "w")
    for line2 in f2:
        line2 = line2.strip()
        if line2.startswith(">"):
            tmp = re.split('\s+',line2,1)
            seqid = tmp[0].split(">")[1]
            if seqid in name:
                data = name[seqid] #改名
                if len(tmp) > 1 and not args.remove_desc:
                    data1 =">"+data+" "+tmp[1]+"\n" #有的fasta文件在序列ID后有描述信息，空格连接
                    f4.write(data1)
                else:
                    data1 =">"+data+"\n"
                    f4.write(data1)
            else:
                f4.write(line2+"\n")
        else:
            f4.write(line2+"\n")
    f2.close()
    f4.close()
                
if args.agp:
#对agp文件的每一行，比一下第一列在不在字典中
    f3 = open(args.agp, "r")
    f5 = open(dout+"/"+args.prefix+".apg", "w")
    for line3 in f3:
        list3 = line3.strip().split("\t")    
        if list3[0] in name:
            data = list3.copy()
            data[0] = name[data[0]]#把第一列改成字典里对应key的value
            data1 = "\t".join(data)+"\n"
            f5.write(data1)
        else:
            f5.write(line3)
    f3.close()
    f5.close()

f1.close()
