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



use constant HELPMSA =><<EOM;


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
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              action = "store", help="Input File[default= \%default]", metavar="character"),
    make_option(c("--title"), type="character", default="bar", 
              action = "store", help="Title[default= \%default]", metavar="character"),
    make_option(c("-b", "--bilv"), type="character", default="200", 
              action = "store", help="ratio[default= \%default]", metavar="character"),
    make_option(c("--outtable"), type="character", default=NULL, 
              action = "store", help="Output Table name[default= \%default]", metavar="character"),
    make_option(c("--outpdf"), type="character", default=NULL, 
              action = "store", help="Output PDF name[default= \%default]", metavar="character")
  );

opt = parse_args(OptionParser(option_list=option_list, usage = "This Script is a test for arguments!"))

if (is.null(opt\$f)){
	print_help(opt_parser)
	stop("Error: invalid input file", call.=FALSE)
}
if (is.null(opt\$outtable)){opt\$outtable=paste(opt\$f,'.xls',sep='')}
if (is.null(opt\$outpdf)){opt\$outpdf=paste(opt\$f,'.pdf',sep='')}



EOM

print "#!/usr/bin/env Rscript\n";
print AUTHORINFO;
print HELPMSA;
