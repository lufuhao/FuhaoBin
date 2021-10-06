#!/usr/bin/env perl

use warnings;
use strict;
###Version: 20212306
use constant AUTHORINFO =><<MEMEMEME;

################# AUTHORS #########################
#  Fu-Hao Lu
#  Post-Doctoral Scientist in Micheal Bevan laboratory
#  Cell and Developmental Department, John Innes Centre
#  Norwich NR4 7UH, United Kingdom
#  E-mail: Fu-Hao.Lu\@jic.ac.uk

MEMEMEME

#argparse.ArgumentParser（）方法参数须知：一般我们只选择用description
#prog=None     - 程序名, (default: sys.argv[0])程序的名字，一般不需要修改，另外，如果你需要在help中使用到程序的名字，可以使用%(prog)s
#description=None,    - help时显示的开始文字
#epilog=None,     - help时显示的结尾文字
#parents=[],        -类型是list，如果这个parser的一些选项跟其他某些parser的选项一样，可以用parents来实现继承，例如parents=[parent_parser]
#formatter_class=argparse.HelpFormatter,     - 自定义帮助信息的格式；三个允许的值： # class argparse.RawDescriptionHelpFormatter 直接输出description和epilog的原始形式（不进行自动换行和消除空白的操作） # class argparse.RawTextHelpFormatter 直接输出description和epilog以及add_argument中的help字符串的原始形式（不进行自动换行和消除空白的操作） # class argparse.ArgumentDefaultsHelpFormatter 在每个选项的帮助信息后面输出他们对应的缺省值，如果有设置的话。这个最常用吧！
#prefix_chars='-',    - 参数的前缀，默认是‘-’，例如-f/–file。有些程序可能希望支持/f这样的选项，可以使用prefix_chars=”/”
#fromfile_prefix_chars=None,     - (default: None)如果你希望命令行参数可以从文件中读取，就可能用到。例如，如果fromfile_prefix_chars=’@’,命令行参数中有一个为”@args.txt”，args.txt的内容会作为命令行参数
#argument_default=None,    - 设置一个全局的选项缺省值，一般每个选项单独设置 (default: None)，所以这个参数用得少，不细说
#conflict_handler='error',     - 不建议使用。这个在极端情况下才会用到，主要是定义两个add_argument中添加的选项的名字发生冲突时怎么处理，默认处理是抛出异常。
#add_help=True    - 是否增加-h/-help选项 (default: True)，一般help信息都是必须的，所以不用设置啦。
#usage - (default: generated)如果你需要修改usage的信息（usage: PROG [-h] [–foo [FOO]] bar [bar …]），那么可以修改这个，一般不要修改。


#add_argument()方法参数须知：
#name or flags...    - 必选，指定参数的形式，一般写两个，一个短参数，一个长参数，看下面的例子”-f”, “–file”
#action	表示值赋予键的方式，这里用到的是bool类型，action意思是当读取的参数中出现指定参数的时候的行为。The basic type of action to be taken when this argument is encountered at the command line. action参数表示值赋予键的方式，这里用到的是bool类型；如果是’count’表示将–verbose标签出现的次数作为verbose的值；'append’表示将每次出现的该便签后的值都存入同一个数组再赋值。
#        argparse内置6种动作可以在解析到一个参数时进行触发：
#           store 保存参数值，可能会先将参数值转换成另一个数据类型。若没有显式指定动作，则默认为该动作。
#           store_const 保存一个被定义为参数规格一部分的值，而不是一个来自参数解析而来的值。这通常用于实现非布尔值的命令行标记。
#           store_ture/store_false 保存相应的布尔值。这两个动作被用于实现布尔开关。
#           append 将值保存到一个列表中。若参数重复出现，则保存多个值。
#           append_const 将一个定义在参数规格中的值保存到一个列表中。
#           version 打印关于程序的版本信息，然后退出


#help		可以写帮助信息 
#required    - 必需参数，通常-f这样的选项是可选的，但是如果required=True那么就是必须的了 #parser.add_argument(“-z”, choices=[‘a’, ‘b’, ‘d’], required=True)
#type   - 指定参数类型。如果希望传进来的参数是指定的类型（例如 float, int or file等可以从字符串转化过来的类型），可以使用 #parser.add_argument(“-x”, type=int) 。
#        这个type类型还可以表示文件操作的类型从而直接进行文件的读写操作。
#       parser.add_argument('file', type=argparser.FileType('r'))    # 读取文件
#       args = parser.parse_args()
#       for line in args.file:
#           print line.strip()
#choices    - 设置参数的范围，如果choice中的类型不是字符串，要指定type, 表示该参数能接受的值只能来自某几个值候选值中，除此之外会报错。choices - 设置参数值的范围，如果choices中的类型不是字符串，记得指定type哦 #parser.add_argument(“-y”, choices=[‘a’, ‘b’, ‘d’])
#nargs    - 指定这个参数后面的value有多少个，默认为1。例如，我们希望使用-n 1 2 3 4，来设置n的值为[1, 2, 3, 4] #parser.add_argument(“-n”, “–num”, nargs=”+”, type=int) # 这里nargs=”+”表示，如果你指定了-n选项，那么-n后面至少要跟一个参数，+表示至少一个,?表示一个或0个,0个或多个 。分析：nargs还可以'*'用来表示如果有该位置参数输入的话，之后所有的输入都将作为该位置参数的值；‘+’表示读取至少1个该位置参数。'?'表示该位置参数要么没有，要么就只要一个。（PS：跟正则表达式的符号用途一致。）
#dest   - 设置这个选项的value解析出来后放到哪个属性中
#default - 如果命令行没有出现这个选项，那么使用default指定的默认值。 #parser.add_argument(“+g”, “++gold”, help=”test test test”,default=”test_gold”)#需要prefix_chars包含”+” 。
#metavar - 参数的名字，在显示 帮助信息时才用到. # parser.add_argument(“-o”, metavar=”OOOOOO”)
#const - A constant value required by some action and nargs selections.

use constant HELPOPT =><<EOP;
# -*- coding: utf-8 -*-
################# Requirements #########################
import sys
import argparse

parser = argparse.ArgumentParser(description = 'this is a description')

parser.add_argument('-s', action='store', required = True, choices = ['test1', 'test2'], dest='simple_value',
        help='Store a simple value')

parser.add_argument('-c', action='store_const', dest='constant_value',
        const='value-to-store',
        help='Store a constant value')

parser.add_argument('-t', action='store_true', default=False,
        dest='boolean_switch',
        help='Set a switch to true')
parser.add_argument('-f', action='store_false', default=False,
        dest='boolean_switch',
        help='Set a switch to false')

parser.add_argument('-a', action='append', dest='collection',
        default=[],
        help='Add repeated values to a list')

parser.add_argument('-A', action='append_const', dest='const_collection',
        const='value-1-to-append',
        default=[],
        help='Add different values to list')
parser.add_argument('-B', action='append_const', dest='const_collection',
        const='value-2-to-append',
        help='Add different values to list')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')

parser.add_argument('--verbose', '-v', action='store_true', help='verbose mode')

results = parser.parse_args()
print 'simple_value     =', results.simple_value
print 'constant_value   =', results.constant_value
print 'boolean_switch   =', results.boolean_switch
print 'collection       =', results.collection
print 'const_collection =', results.const_collection

                                                                                   # 参数字典的verbose建对应的值为True，而help参数用于描述--verbose参数的用途或意义。
# 将变量以标签-值的字典形式存入args字典
if args.verbose:
    print "Verbose mode on!"
else:
    print "Verbose mode off!"
    

EOP


print "#!/usr/bin/env python3\n";

print HELPOPT;
print AUTHORINFO;
