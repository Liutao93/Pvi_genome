############################################################################
#北京组学生物科技有限公司
#author whj
#date 2022.10.27
#version 1.0
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
parser <- ArgumentParser(description='ConsensusCluster analysis')
parser$add_argument( "-i", "--input_divsum", type="character",required=F,
                     help="input the immune martix [required]",metavar="filepath")
parser$add_argument( "-g", "--genomesize", type="integer",required=T,
                     help="input the Output file prefix",metavar="genomesize")
parser$add_argument( "-p", "--divsumPrefix", type="character",required=F,default="prefixname",
                     help="input the Output file prefix",metavar="outdir")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default %(default)s]",
                     metavar="outdir")
parser$add_argument( "-H", "--height", type="double", default=4.5,
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
#divsyum data process
#————————————————————————————————————————————————————————————————————
TE_divsum1 = read.table(opt$input_divsum,header = T,sep = " ",check.names = F) %>%
  pivot_longer(-Div, names_to = "Classification",values_to = "size") %>%
  mutate(percent = (size/opt$genomesize) * 100)
head(TE_divsum1)

TE_divsum1$class=str_split(TE_divsum1$Classification,"/",simplify = T)[,1]
head(TE_divsum1)
TE_divsum1[which(TE_divsum1$class == 'RC'|TE_divsum1$class=='RC.'),'class'] <- 'DNA'
TE_divsum1=TE_divsum1[which(TE_divsum1$class == 'DNA'|TE_divsum1$class == 'LTR'|TE_divsum1$class == 'LINE'|TE_divsum1$class == 'SINE'),]
TE_divsum1$class=factor(TE_divsum1$class, levels=c('LTR', 'DNA', 'LINE', 'SINE'))
TE_divsum1[which(TE_divsum1$Div >=50),'Div']=50
#columnar graph
gra_2<-ggplot(data = TE_divsum1, aes(x = Div, y = percent)) +
  geom_col(aes(fill = class), width = 0.8, color = F) +
  scale_x_continuous(breaks = c(0,10,20,30,40,50),
                     labels = c(0,10,20,30,40,'>50')) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Kimura subsitution levels (CpG adjusted)",
       y = "Percent of the genome(%)", fill = NULL) +
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.position = c(0.95,0.95),
    legend.justification = c("right", "top"))
ggsave(gra_2,filename = paste0(opt$outdir,'/',opt$divsumPrefix,'.png'),width = opt$width,height = opt$height)
ggsave(gra_2,filename = paste0(opt$outdir,'/',opt$divsumPrefix,'.pdf'),width = opt$width,height = opt$height)
