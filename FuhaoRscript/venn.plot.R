#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

InData <- args[1]
x<- read.table(InData, header=TRUE, sep="\t")
#library(venn)
#library(stringr)
#setEPS()
#postscript(paste(args[1], ".eps", sep=""), width=8, height=8, pointsize=10)

#venn(list(D10G=x[x[2]==1,1],D27G=x[x[3]==1,1], Leaf=x[x[4]==1,1], Root=x[x[5]==1,1], Seedling=x[x[6]==1,1]), ilabels=TRUE, ellipse=TRUE, zcolor="style", size=8, borders = TRUE)
#venn(list(DEG=x[x[2]==1,1],dCDS=x[x[3]==1,1], dUTR3=x[x[4]==1,1], dUTR5=x[x[5]==1,1], dATAC=x[x[6]==1,1]), ilabels=TRUE, zcolor="style", borders = TRUE, main="ATAC-seq vs DEGs")

#dev.off()

#tiff(filename = paste(args[1], ".tif", sep=""),
#     width = 8, height = 8, units = "cm", pointsize = 10,
#     compression = "lzw", bg = "white", res = 600, family = "Arial")
#venn(list(DEG=x[x[2]==1,1],dCDS=x[x[3]==1,1], dUTR3=x[x[4]==1,1], dUTR5=x[x[5]==1,1], dATAC=x[x[6]==1,1]), ilabels=TRUE, zcolor="style", borders = TRUE, main="ATAC-seq vs DEGs")
#dev.off()



library(Vennerable)
#tiff(filename = paste(args[1], ".1.tif", sep=""),
#     width = 8, height = 8, units = "cm", pointsize = 10,
#     compression = "lzw", bg = "white", res = 600, family = "Arial")
setEPS()
postscript(paste(args[1], ".1.eps", sep=""), width=8, height=8, pointsize=10)
venndata <- Venn(list("DEG"=x[x[2]==1,1], "dATAC"=x[x[6]==1,1]))
plot(venndata, doWeights = TRUE, doEuler=FALSE, type="squares")
dev.off()

#tiff(filename = paste(args[1], ".2.tif", sep=""),
#     width = 8, height = 8, units = "cm", pointsize = 10,
#     compression = "lzw", bg = "white", res = 600, family = "Arial")
setEPS()
postscript(paste(args[1], ".2.eps", sep=""), width=8, height=8, pointsize=10)
venndata <- Venn(list("dUTR5"=x[x[5]==1,1], "dUTR3"=x[x[4]==1,1], "dCDS"=x[x[3]==1,1]))
plot(venndata, doWeights = TRUE, doEuler=TRUE, type="circles")
dev.off()
