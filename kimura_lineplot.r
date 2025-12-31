############################################################################
#北京组学生物科技有限公司
#author whj
#date 2022.10.27
#version 1.0

#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.h5.xeknow.com/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.h5.xeknow.com/s/26edpc

###############################################################################################

library("tidyverse")
library("ggplot2")
library('argparse')
p="argparse"
parser <- ArgumentParser(description='plot analysis')
parser$add_argument( "-i", "--input_divsum", type="character",required=T,
                     help="input the immune martix [required]",metavar="filepath")
parser$add_argument( "-p", "--outputPrefix", type="character",required=T,default="line",
                     help="input the Output file prefix,default=%(default)s")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default %(default)s]",
                     metavar="outdir")
parser$add_argument( "-H", "--height", type="double", default=5,
                     help="the height of pic inches  [default %(default)s]",
                     metavar="height")
parser$add_argument("-W", "--width", type="double", default=5,
                    help="the width of pic inches [default %(default)s]",
                    metavar="width")
opt <- parser$parse_args()
if( !file.exists(opt$outdir) ){
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
}

#————————————————————————————————————————————————————————————————————
# gff3 data process
#————————————————————————————————————————————————————————————————————
allTE_RM = read.table(opt$input_divsum,skip = 7,header = F)
# allTE_RM = read.table('/share/nas5/wanghj/data/work_2/genome.identity',skip = 7,header = F)
allTE_RM$class=  str_split(allTE_RM$V1, "/", simplify = T)[,1]
head(allTE_RM,20)

allTE_RM[which(allTE_RM$class == 'RC'|allTE_RM$class=='RC?'),'class'] <- 'DNA'
allTE_RM=allTE_RM[which(allTE_RM$class == 'DNA'|allTE_RM$class == 'LTR'|
                          allTE_RM$class == 'LINE'|allTE_RM$class == 'SINE'),]
allTE_RM$class=factor(allTE_RM$class, levels=c('LTR', 'DNA', 'LINE', 'SINE'))
# allTE_RM=allTE_RM[which(allTE_RM$Identity>=0),]

#allTE_RM$V5
allTE_RM$Identity=as.double(allTE_RM$V5)/100
#line graph
gra_1<-ggplot(data = allTE_RM,aes(x=Identity))+
  scale_color_brewer(palette = "Set1")+
  # geom_density(stat='bin',bins = 30,aes(color = class)) +
  geom_density(aes(x=Identity,y=..count../45,color = class))+
  scale_x_continuous(expand = c(0, 0),limits = c(-0.1, 1.1)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 730)) +
  labs(x = "Weighted average Kimura divergence",y='count') +
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position = c(0.9,0.8))
ggsave(gra_1,filename = paste0(opt$outdir,'/',opt$outputPrefix,'.png'),width = opt$width,height = opt$height)
ggsave(gra_1,filename = paste0(opt$outdir,'/',opt$outputPrefix,'.pdf'),width = opt$width,height = opt$height)
