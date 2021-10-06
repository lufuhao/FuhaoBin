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

#argparse.ArgumentParser��������������֪��һ������ֻѡ����description
#prog=None     - ������, (default: sys.argv[0])��������֣�һ�㲻��Ҫ�޸ģ����⣬�������Ҫ��help��ʹ�õ���������֣�����ʹ��%(prog)s
#description=None,    - helpʱ��ʾ�Ŀ�ʼ����
#epilog=None,     - helpʱ��ʾ�Ľ�β����
#parents=[],        -������list��������parser��һЩѡ�������ĳЩparser��ѡ��һ����������parents��ʵ�ּ̳У�����parents=[parent_parser]
#formatter_class=argparse.HelpFormatter,     - �Զ��������Ϣ�ĸ�ʽ�����������ֵ�� # class argparse.RawDescriptionHelpFormatter ֱ�����description��epilog��ԭʼ��ʽ���������Զ����к������հ׵Ĳ����� # class argparse.RawTextHelpFormatter ֱ�����description��epilog�Լ�add_argument�е�help�ַ�����ԭʼ��ʽ���������Զ����к������հ׵Ĳ����� # class argparse.ArgumentDefaultsHelpFormatter ��ÿ��ѡ��İ�����Ϣ����������Ƕ�Ӧ��ȱʡֵ����������õĻ��������ðɣ�
#prefix_chars='-',    - ������ǰ׺��Ĭ���ǡ�-��������-f/�Cfile����Щ�������ϣ��֧��/f������ѡ�����ʹ��prefix_chars=��/��
#fromfile_prefix_chars=None,     - (default: None)�����ϣ�������в������Դ��ļ��ж�ȡ���Ϳ����õ������磬���fromfile_prefix_chars=��@��,�����в�������һ��Ϊ��@args.txt����args.txt�����ݻ���Ϊ�����в���
#argument_default=None,    - ����һ��ȫ�ֵ�ѡ��ȱʡֵ��һ��ÿ��ѡ������� (default: None)��������������õ��٣���ϸ˵
#conflict_handler='error',     - ������ʹ�á�����ڼ�������²Ż��õ�����Ҫ�Ƕ�������add_argument����ӵ�ѡ������ַ�����ͻʱ��ô����Ĭ�ϴ������׳��쳣��
#add_help=True    - �Ƿ�����-h/-helpѡ�� (default: True)��һ��help��Ϣ���Ǳ���ģ����Բ�����������
#usage - (default: generated)�������Ҫ�޸�usage����Ϣ��usage: PROG [-h] [�Cfoo [FOO]] bar [bar ��]������ô�����޸������һ�㲻Ҫ�޸ġ�


#add_argument()����������֪��
#name or flags...    - ��ѡ��ָ����������ʽ��һ��д������һ���̲�����һ��������������������ӡ�-f��, ���Cfile��
#action	��ʾֵ������ķ�ʽ�������õ�����bool���ͣ�action��˼�ǵ���ȡ�Ĳ����г���ָ��������ʱ�����Ϊ��The basic type of action to be taken when this argument is encountered at the command line. action������ʾֵ������ķ�ʽ�������õ�����bool���ͣ�����ǡ�count����ʾ���Cverbose��ǩ���ֵĴ�����Ϊverbose��ֵ��'append����ʾ��ÿ�γ��ֵĸñ�ǩ���ֵ������ͬһ�������ٸ�ֵ��
#        argparse����6�ֶ��������ڽ�����һ������ʱ���д�����
#           store �������ֵ�����ܻ��Ƚ�����ֵת������һ���������͡���û����ʽָ����������Ĭ��Ϊ�ö�����
#           store_const ����һ��������Ϊ�������һ���ֵ�ֵ��������һ�����Բ�������������ֵ����ͨ������ʵ�ַǲ���ֵ�������б�ǡ�
#           store_ture/store_false ������Ӧ�Ĳ���ֵ������������������ʵ�ֲ������ء�
#           append ��ֵ���浽һ���б��С��������ظ����֣��򱣴���ֵ��
#           append_const ��һ�������ڲ�������е�ֵ���浽һ���б��С�
#           version ��ӡ���ڳ���İ汾��Ϣ��Ȼ���˳�


#help		����д������Ϣ 
#required    - ���������ͨ��-f������ѡ���ǿ�ѡ�ģ��������required=True��ô���Ǳ������ #parser.add_argument(��-z��, choices=[��a��, ��b��, ��d��], required=True)
#type   - ָ���������͡����ϣ���������Ĳ�����ָ�������ͣ����� float, int or file�ȿ��Դ��ַ���ת�����������ͣ�������ʹ�� #parser.add_argument(��-x��, type=int) ��
#        ���type���ͻ����Ա�ʾ�ļ����������ʹӶ�ֱ�ӽ����ļ��Ķ�д������
#       parser.add_argument('file', type=argparser.FileType('r'))    # ��ȡ�ļ�
#       args = parser.parse_args()
#       for line in args.file:
#           print line.strip()
#choices    - ���ò����ķ�Χ�����choice�е����Ͳ����ַ�����Ҫָ��type, ��ʾ�ò����ܽ��ܵ�ֵֻ������ĳ����ֵ��ѡֵ�У�����֮��ᱨ��choices - ���ò���ֵ�ķ�Χ�����choices�е����Ͳ����ַ������ǵ�ָ��typeŶ #parser.add_argument(��-y��, choices=[��a��, ��b��, ��d��])
#nargs    - ָ��������������value�ж��ٸ���Ĭ��Ϊ1�����磬����ϣ��ʹ��-n 1 2 3 4��������n��ֵΪ[1, 2, 3, 4] #parser.add_argument(��-n��, ���Cnum��, nargs=��+��, type=int) # ����nargs=��+����ʾ�������ָ����-nѡ���ô-n��������Ҫ��һ��������+��ʾ����һ��,?��ʾһ����0��,0������ ��������nargs������'*'������ʾ����и�λ�ò�������Ļ���֮�����е����붼����Ϊ��λ�ò�����ֵ����+����ʾ��ȡ����1����λ�ò�����'?'��ʾ��λ�ò���Ҫôû�У�Ҫô��ֻҪһ������PS����������ʽ�ķ�����;һ�¡���
#dest   - �������ѡ���value����������ŵ��ĸ�������
#default - ���������û�г������ѡ���ôʹ��defaultָ����Ĭ��ֵ�� #parser.add_argument(��+g��, ��++gold��, help=��test test test��,default=��test_gold��)#��Ҫprefix_chars������+�� ��
#metavar - ���������֣�����ʾ ������Ϣʱ���õ�. # parser.add_argument(��-o��, metavar=��OOOOOO��)
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

                                                                                   # �����ֵ��verbose����Ӧ��ֵΪTrue����help������������--verbose��������;�����塣
# �������Ա�ǩ-ֵ���ֵ���ʽ����args�ֵ�
if args.verbose:
    print "Verbose mode on!"
else:
    print "Verbose mode off!"
    

EOP


print "#!/usr/bin/env python3\n";

print HELPOPT;
print AUTHORINFO;
