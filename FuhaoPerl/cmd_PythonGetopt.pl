#!/usr/bin/env perl

use warnings;
use strict;
###Version: 20200803
use constant AUTHORINFO =><<MEMEMEME;

################# AUTHORS #########################
#  Fu-Hao Lu
#  Post-Doctoral Scientist in Micheal Bevan laboratory
#  Cell and Developmental Department, John Innes Centre
#  Norwich NR4 7UH, United Kingdom
#  E-mail: Fu-Hao.Lu\@jic.ac.uk

MEMEMEME



use constant HELPOPT =><<EOP;
# -*- coding: utf-8 -*-
################# Requirements #########################
import getopt
import sys

def ArgsParser(argv):
  functionname = 'ArgsParse'
  partdate = '20200803'
   
  try:
    opts, args = getopt.getopt(argv, "hf:o:", ["help", "fasta=", "output="]) 
　　　　　#表示参数选项有：-h, -f, -o, --help, --fasta, --output，它们相互对应；该方法的返回值有两个元素: 第一个是(opt, value)元组的列表，第二个是一般参数列表，包含那些没有 '-' 或 '--' 的参数
  except getopt.GetoptError:
    print('Error: ***.py -f <fasta file> -o <output>')
    print('  or: ***.py --fasta=<fasta file> --output=<partdate>')
    sys.exit(2)
   
  for opt, arg in opts:  #依次获取列表中的元组项
    if opt in ("-h", "--help"):
      print('***.py -f <fasta file> -o <output>')
      print('or: ***.py --fasta=<functionname> --output=<partdate>')
      sys.exit()
    elif opt in ("-i", "--fasta"):
      fasta_file = arg
    elif opt in ("-o", "--output"):
      output = arg
  print('-----------------------------------------------------------------------')
  print(opts) #元组构成的列表
  print(args) #args指的是不用 '-'或 '--'传递的参数，这里没有传递，所以为空
  print('Input Fasta: ', fasta_file)
  print('Ouput', output)
 
if __name__ == '__main__':
  ArgsParser(sys.argv[1:]) #因为sys.argv[0]是脚本名称



EOP


print "#!/usr/bin/env python3\n";
print AUTHORINFO;
print HELPOPT;
