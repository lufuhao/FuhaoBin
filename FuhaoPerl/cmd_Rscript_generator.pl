#!/usr/bin/env perl
use warnings;
use strict;
###Version: 20200801
use constant AUTHORINFO =><<MEMEMEME;

################# AUTHORS #########################
#  Fu-Hao Lu
#  Post-Doctoral Scientist in Micheal Bevan laboratory
#  Cell and Developmental Department, John Innes Centre
#  Norwich NR4 7UH, United Kingdom
#  E-mail: Fu-Hao.Lu\@jic.ac.uk

MEMEMEME



use constant HELPSIM =><<EOH;

args <- commandArgs(trailingOnly = TRUE)
filename  =args[1]
### can not define default, or omit parameters
commandArgs(T)

EOH



use constant HELPOPT =><<EOP;

################# Requirements #########################
# getopt
#     install.packages("getopt")
#     useage: getopt(spec = NULL, opt = commandArgs(TRUE),command = get_Rscript_filename(), usage = FALSE,debug = FALSE)
#         spec: a matrix
#             col1: parameter name	"first_arg"
#             col2: parameter short name	"f"
#             col3: parameter pattern	0 no param 1 both 2 must para 
#             col4: parameter class	logical; integer; double; complex; character; numeric
#             col5: note	"this arg is the first argument"
#         opt: parameter source, commandArg() by default
#         usage: help message, FALSE by default


library(getopt)

# 首先将第一个参数的具体信息进行描述
spec <- matrix(
# 每行五个，第五个可选，也就是说第五列可以不写
# byrow 按行填充矩阵的元素
# ncol  每行填充五个元素
  c("first",  "f", 2, "integer", "This is first!",
    "second", "s", 1, "character",  "This is second!",
    "third",  "t", 2, "double",  "This is third!",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)

# 使用getopt方法，注意这里的usage默认是关闭的，这里不能打开
opt <- getopt(spec=spec)

# 这个时候需要将usage参数打开，这样getpot()就会返回一个特定写法的帮助文件
# 当然你也可以自己写帮助，然后将if判断语句中的打印的信息换成你自己定义的帮助信息
if( !is.null(opt\$help) || is.null(opt\$first) || is.null(opt\$third) ){
    # ... 这里你也可以自定义一些东放在里面
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

# opt实际上就是一个列表，直接使用\$来索引到对应的参数的值
print(opt\$first)
print(opt\$second)
print(opt\$third)



EOP



use constant HELPMAS1 =><<EOM;


is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 


################# Requirements #########################
# optparse
#     install.packages("optparse")
#         OptionParser(usage = "usage: \%prog [options]", option_list = list(),
#                      add_help_option = TRUE, prog = NULL, description = "",
#                      epilogue = "")
packages.required <- c("optparse", "dplyr", "ggplot2", "grid", "RColorBrewer", "SuppDists", "scales", "tidyr")



if (!requireNamespace("BiocManager", quietly = T))
        install.packages("BiocManager")
invisible(library("BiocManager", character.only = T, quietly = T, warn.conflicts = FALSE))
packages.absense <- setdiff(packages.required, rownames(installed.packages()))
if (length(packages.absense) > 0) {
	install.packages(packages.absense)
}
#lapply(packages.required, library, quiet=T)
#lapply(packages.required, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))
invisible(lapply(packages.required, function(x) require(x, character.only = T, quietly = T, warn.conflicts = FALSE)))

if(! require("optparse")) install.packages("optparse")
library(optparse)
if(! require("LinkageMapView")) install.packages("LinkageMapView")
library(LinkageMapView)
if(! require("ComplexHeatmap"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)



packages <- c("ggplot2", "dplyr", "lavaan", "plyr", "cowplot", "rmarkdown", "readr", "caTools", "bitops")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	install.packages(setdiff(packages, rownames(installed.packages())))
}



### Alternative
for (pkg in cran_packages){
	if (! require(pkg,character.only=T) ) {
		install.packages(pkg,ask = F,update = F)
		require(pkg,character.only=T) 
	}
}
# first prepare BioManager on CRAN
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
# use BiocManager to install
for (pkg in Biocductor_packages){
	if (! require(pkg,character.only=T) ) {
		BiocManager::install(pkg,ask = F,update = F)
		require(pkg,character.only=T) 
	}
}
for (pkg in c(Biocductor_packages,cran_packages)){
	require(pkg,character.only=T) 
}
if(!require(tinyarray))devtools::install_local("tinyarray-master.zip",upgrade = F)






#注意，这个模块不用加上-h的flag，不然会报错
option_list = list(
#    make_option(c("-f", "--file"), type="logical/integer/double/complex/character", default=NULL, 
#              action = "store", dest=file, help="Input File[default= \%default]", metavar="character"),
    make_option(c("-f", "--file"), type="character", default=NULL, 
              action = "store", dest=file, help="Input File[default= \%default]", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="MyOut", 
              action = "store", dest="out", help="Output File prefix [default= \%default]", metavar="PREFIX"),
    make_option(c("--eps"), type="logical", default=FALSE, 
              action = "store_true", help="Output EPS [default= \%default]", metavar=NULL),
    make_option(c("--svg"), type="logical", default=FALSE, 
              action = "store_true", help="Output SVG [default= \%default]", metavar=NULL),
    make_option(c("--tif"), type="logical", default=FALSE, 
              action = "store_true", help="Output TIFF [default= \%default]", metavar=NULL),
    make_option(c("--pdf"), type="logical", default=FALSE, 
              action = "store_true", help="Output TIFF [default= \%default]", metavar=NULL),
    make_option(c("--width"), type="double", default=8, 
              action = "store", help="Output TIFF/SVG/EPS/PDF width (cm) [default= \%default]", metavar="float"),
    make_option(c("--height"), type="double", default=8, 
              action = "store", help="Output TIFF/SVG/EPS/PDF height (cm) [default= \%default]", metavar="float"),
    make_option(c("--resolution"), type="double", default=600, 
              action = "store", help="Output TIFF resolution [default= \%default]", metavar="float"),
    make_option(c("--pointsize"), type="double", default=10, 
              action = "store", help="Output TIFF/SVG/EPS/PDF [default= \%default]", metavar="float"),
    make_option(c("--family"), type="character", default="ArialMT", 
              action = "store", help="Output TIFF/SVG/EPS/PDF font family [default= \%default]", metavar="character"),
    make_option(c("--verbose"), type="logical", default=FALSE, 
              action = "store_true", help="Print detailed info [default= \%default]", metavar=NULL)
  );

#opt = parse_args(OptionParser(option_list=option_list, usage = "This Script is a test for arguments!"))
opt_parser=OptionParser(option_list=option_list, usage = "This Script is  to display blast results!")
opt = parse_args(opt_parser)

if (is.null(opt\$file)){
	print_help(opt_parser)
	stop("Error: invalid input file", call.=FALSE)
}
if (is.null(opt\$outtable)){opt\$outtable=paste(opt\$file,'.xls',sep='')}
if (is.null(opt\$outpdf)){opt\$outpdf=paste(opt\$file,'.pdf',sep='')}


################# input and output #########################
print (paste0("Input File: ", opt\$file))
print (paste0("Output File: ", opt\$out))


### BitMap
#bitmap(file = "test1.jpeg", type = "jpeg", res = 1200) plot(1:22, pch = 1:22, cex = 1:3, col = 1:5)
#dev.off()
#bmp(filename = "Rplot%03d.bmp", width = 480, height = 480, units = c("px", "cm"), pointsize = 12, bg = "white", res = NA, ..., family = "", restoreConsole = TRUE, type = c("cairo", "Xlib", "quartz"), antialias)

### EPS
#setEPS()
#postscript(OutEps, width=8, height=2, pointsize=10)
#postscript(file="heatmaps_in_r.eps", onefile=FALSE, horizontal=FALSE, width = 1000, height = 2000)
#dev.off()

### JPEG
#jpeg(filename = "Rplot%03d.jpeg", width = 480, height = 480, units = "px", pointsize = 12, quality = 75, bg = "white", res = NA, ..., family = "", restoreConsole = TRUE, type = c("cairo", "Xlib", "quartz"), antialias)
#dev.off()

### PDF
#pdf(file="myplot.pdf")
#dev.off()
###PDFs are 7x7 inches by default, and each new plot is on a new page. The size can be changed:
### 6x3 inches
#pdf("plots.pdf", width=6, height=3)
### 10x6 cm
pdf("plots.pdf", width=10/2.54, height=6/2.54)
###If you want to edit your file in a vector editor like Inkscape or Illustrator, some of the plotting point objects might look like letters instead of circles, squares, etc. To avoid this problem:
#pdf("plots.pdf", useDingbats=FALSE)

### PNG
### By default, the graphs are 480x480 pixels in size, at a resolution of 72 dpi (6.66x6.66 inches).
#png(file="myplot.png", bg="transparent", width=480, height=240, res=120)
#dev.off()
#png(filename = "Rplot%03d.png", width = 480, height = 480, units = "px", pointsize = 12, bg = "white", res = NA, ..., family = "", restoreConsole = TRUE, type = c("cairo", "cairo-png", "Xlib", "quartz"), antialias)

### PS
#postscript("a.ps")
#dev.off()


### SVG
#svg(filename="Std_SVG.svg", width=5, height=4, pointsize=12)
#dev.off()

### TIFF
#tiff(file = "C:/test1.tiff", res = 300, width = 2400, height = 2400, compression = "lzw")
#plot(1:22, pch = 1:22, cex = 1:3, col = 1:5)
#dev.off()
#tiff(filename = "Rplot%03d.tiff", width = 480, height = 480, units = "px", pointsize = 12, compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"), bg = "white", res = NA, ..., family = "Arial", restoreConsole = TRUE, type = c("cairo", "Xlib", "quartz"), antialias)
#tiff(filename = paste(args[1], ".tif", sep=""), width = 8, height = 8, units = "cm", pointsize = 12, compression = "lzw", bg = "white", res = 600, family = "Arial")

if (opt\$svg) {
	svg(filename=paste(fig.out.pfx, ".svg", sep=""), width=opt\$width/cm(1), height=opt\$height/cm(1), pointsize=opt\$pointsize, bg = "white", family=opt\$family)
	My.Plot
	dev.off()
}
if (opt\$pdf) {
	pdf(file=paste(fig.out.pfx, ".pdf", sep=""), width=opt\$width/cm(1), height=opt\$height/cm(1), pointsize=opt\$pointsize, family=opt\$family, bg = "white")
	My.Plot
	dev.off()
}
if (opt\$tif) {
	tiff(filename = paste(fig.out.pfx, ".tif", sep=""), width=opt\$width, height=opt\$height, units = "cm", pointsize=opt\$pointsize, compression = "lzw", bg = "white", res=opt\$resolution, family=opt\$family)
	My.Plot
	dev.off()
}
if (opt\$eps) {
	setEPS()
	postscript(paste(fig.out.pfx, ".eps", sep=""), width=opt\$width/cm(1), height=opt\$height/cm(1), pointsize=opt\$pointsize, family=opt\$family, bg = "white")
	My.Plot
	dev.off()
}

EOM

use constant HELPMAS2 =><<EOM2;


is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 


################# Requirements #########################
# optparse
#     install.packages("optparse")
#         OptionParser(usage = "usage: \%prog [options]", option_list = list(),
#                      add_help_option = TRUE, prog = NULL, description = "",
#                      epilogue = "")
if(! require("optparse")) install.packages("optparse")
library(optparse)
if(! require("LinkageMapView")) install.packages("LinkageMapView")
library(LinkageMapView)


#注意，这个模块不用加上-h的flag，不然会报错
library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true", 
                default=TRUE, help="Print extra output [default]")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false", 
                    dest="verbose", help="Print little output")
parser <- add_option(parser, c("-c", "--count"), type="integer", default=5, 
                help="Number of random normals to generate [default \%default]",
                metavar="number")
parse_args(parser, args = c("--quietly", "--count=15"))

## \$help
## [1] FALSE
## 
## \$verbose
## [1] FALSE
## 
## \$count
## [1] 15



EOM2

print "#!/usr/bin/env Rscript\n";
print AUTHORINFO;
print HELPMAS1;
