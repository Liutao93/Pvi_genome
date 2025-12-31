############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2022.12.01
#version 1.0

#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.h5.xeknow.com/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.h5.xeknow.com/s/26edpc

###############################################################################################

library('argparse')
p="argparse"
parser <- ArgumentParser(description='plot  analysis')
parser$add_argument( "-i", "--input", type="character",required=T,
                     help="input the immune martix [required]",metavar="filepath")
parser$add_argument( "-p", "--prefix", type="character",required=F,default="line",
                     help="input the Output file prefix,default=%(default)s")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default %(default)s]",
                     metavar="outdir")
parser$add_argument( "-H", "--height", type="double", default=5,
                     help="the height of pic inches  [default %(default)s]",
                     metavar="height")
parser$add_argument("-W", "--width", type="double", default=8,
                    help="the width of pic inches [default %(default)s]",
                    metavar="width")
opt <- parser$parse_args()
if( !file.exists(opt$outdir) ){
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
}



library("tidyverse")
library("ggplot2")
#————————————————————————————————————————————————————————————————————
# gff3 data process
#————————————————————————————————————————————————————————————————————
alldata = read.table(opt$input,header = F,sep="\t",check.names = F,comment.char = "")

colnames(alldata)<-c("type","value","Species")

gene_length=subset(alldata,type=="gene_length" & value<20000)
cds_length=subset(alldata,type=="cds_length" &value<6000)
cds_number=subset(alldata,type=="cds_number" &value<50)
introns_length=subset(alldata,type=="introns_length" &value<5000)
exons_length=subset(alldata,type=="exons_length" &value<600)
exon_number=subset(alldata,type=="exon_number" &value<50)

subsetData=rbind(gene_length,cds_length,cds_number,introns_length,exons_length,exon_number)

subsetData$type=factor(subsetData$type,levels = c("gene_length","cds_length","cds_number","introns_length","exons_length","exon_number"),order=T)


#opt=list()
#opt$input="all_species_gene_structure.tsv"
#line graph
g<-ggplot(data = subsetData,aes(x=value))+
  stat_density(aes(x=value,y=..count.., colour=Species),
               geom="line",position="identity")+
  facet_wrap(vars(type), scales = "free")+
  scale_color_brewer(palette = "Set2")+
  labs(x = "",y='Count') +
  theme_bw()+
  theme(panel.grid=element_blank())
  #theme(legend.position = c(0.9,0.8))
ggsave(g,filename = paste0(opt$outdir,'/',opt$prefix,'.png'),width = opt$width,height = opt$height,dpi=300)
ggsave(g,filename = paste0(opt$outdir,'/',opt$prefix,'.pdf'),width = opt$width,height = opt$height)
