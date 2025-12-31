############################################################################
#北京组学生物科技有限公司
#author whj
#date 2022.12.03
#version 1.0
#学习R课程：
#python 学习课程:https://zzw.xet.tech/s/2bdd89
###############################################################################################

import os, re,argparse,sys

parser = argparse.ArgumentParser(description="help")
parser.add_argument('-i', "--inputfile", required=True,help=" input stat file ",)
parser.add_argument('-o', "--outputdir",default=os.getcwd(), help="specify the output file dir ,default %(default)s",)
parser.add_argument('-p','--prefix',default='genome',help=' specify the output file prefix , default %(default)s',)
#parser.add_argument("--upper-case",default=False,action="store_true",help=" help ", )

args = parser.parse_args()

f = args.inputfile
ip=args.outputdir+'/'+args.prefix+'.gff'
# data_operate.updateFile(0,f,'qweadsfasrqWR')
with open(ip, 'w') as g:
    for i in open(f):
        #i = re.split("\t",i)
        i=i.split()
        if len(i) >= 15:
            if 'Low' in i[10] or 'Simple' in i[10] or "Satellite" in i[10]: continue
            j = (i[4],      #1 seqid
                 'TRF' if 'Low' in i[10] or 'Simple' in i[10] else 'RepeatMasker',       #2 source
                 'TandemRepeat' if 'Low' in i[10] or 'Simple' in i[10] else 'Transposon',       #3 type
                 i[5],      #4 start
                 i[6],      #5 end
                 i[0],      #6 score
                 i[8] if i[8] != 'C' else '-',      #7 strand
                 '.',       #8 phase
                 ';'.join(['ID=TP' + i[14].zfill(7), 'class=' + i[10], 'matching=' + i[9]]))
            g.write('\t'.join(j) + '\n')
