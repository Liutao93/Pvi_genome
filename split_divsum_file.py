
############################################################################
#北京组学生物科技有限公司
#author whj
#date 2022.10.27
#version 1.0
#学习R课程：
#python 学习课程:https://zzw.xet.tech/s/2bdd89
###############################################################################################

import os,re,argparse

parser = argparse.ArgumentParser(description="help")
parser.add_argument('-i', "--inputfile", help="\033[1;36m input stat file \033[0m",)
parser.add_argument('-o', "--outputdir",default=os.getcwd(), help="\033[1;36m specify the output file dir,default is current \033[0m",)
parser.add_argument('-p','--outFilePrefix',default='genome',help='\033[1;36m specify the output file prefix,defaul demo \033[0m',)
parser.add_argument("--upper-case",default=False,action="store_true",help="\033[1;36m help \033[0m", )

args = parser.parse_args()


if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)


def updateFile(fun,file, old_str, *new_str):
    #1 is replace，0 is delect
    with open(file, "r") as f1, open("%s.bak" % file, "w") as f2:
        if fun == 1:
            for i in f1:
                f2.write(re.sub(old_str, *new_str, i))
        if fun == 0:
            for i in f1:
                if old_str not in i and i!='\n':
                    f2.write(i)
    os.remove(file)
    os.rename("%s.bak" % file, file)


arr = open(args.inputfile).read().split('\n\n')

file_d=args.outputdir+'/'+args.outFilePrefix+'.divsum'

with open(file_d, 'w') as f:
    f.write(arr[1])

updateFile(1,file_d,'/','_')
updateFile(0,file_d,'Coverage')

