
#!/usr/bin/Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2021.05.19
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.h5.xeknow.com/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.h5.xeknow.com/s/26edpc


#设置镜像，
local({r <- getOption("repos")  
r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos=r)}) 



#2 安装bioconductor常用包
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")



p="argparse"

if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
  install.packages(p,  warn.conflicts = FALSE)
  suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
}


#######命令行参数设置################################################################################
parser <- ArgumentParser(description='merge tsv files by rownames: https://www.omicsclass.com/article/1494 ')

parser$add_argument( "-i", "--input", type="character",nargs="+",required=T,
                     help="input gct format file [required]",
                     metavar="input")
parser$add_argument( "-b", "--by", type="character",nargs="+",required=F, default="barcode",
                     help="character vector of variable names to join by. If omitted,will match on all common variables. ",
                     metavar="by")
parser$add_argument( "-t", "--type", type="character",required=F, default="left",
                     help="type of join: left (default), right, inner or full ",
                     metavar="type")

parser$add_argument( "-p", "--prefix", type="character", default="demo",
                     help="output file name prefix [default demo]",
                     metavar="prefix")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default cwd]",
                     metavar="outdir")

opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
  }
}

#######################################################################################

package_list <- c("plyr")
# 判断R包加载是否成功如果加载不成功自动安装
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,  warn.conflicts = FALSE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

###########################################################################

filenames=opt$input
#opt
datalist = lapply(filenames, function(x){
                read.table(file=x,header=T,sep = "\t",check.names = F)
                  })

by=NULL
if (length(opt$by)>1){
  for(i in 1:length(opt$by)){
    if(i==1){
      by=opt$by[1]
      next
    }else{
      if(by %in% colnames(datalist[[i]])){
        
      }else{
        datalist[[i]][,by]=datalist[[i]][,opt$by[i]]
      }
      
      
    }
  }
  
  
}else{
  by=opt$by
}


dd=join_all(datalist,by=by, type=opt$type)
write.table(dd,file =paste0(opt$outdir,"/",opt$prefix,".tsv"),sep="\t",quote = F,row.names = F)
