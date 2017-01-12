#!/usr/bin/env python
# -*- coding: utf-8 -*-
#http://techoverflow.net/blog/2013/11/05/gtf2gff.py-a-replacement-for-gtf2gff.pl/
"""
gtf2gff.py -- A script to convert GTF to GFF files.
... and a better replacement for gtf2gff.pl

Version 1.0
"""
# Python 2.7
from __future__ import with_statement
import argparse
import sys

__author__    = "Uli Koehler & Anton Smirnov"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__   = "Apache v2.0"

class GTFException(Exception):
    pass

def gtf2gff(infileName, outfileName, startindex, endindex, program):
    with open(infileName,"r") as infile, open(outfileName,"w") as outfile:
        genId = 0
        for line in infile:
            line = line.strip()
            if not line: continue
            words = line.split("\t")
            if len(words) != 9:
                raise GTFException("Encountered %d columns instead of the expected 9 in line: '%s'" % (len(words), line))
            if words[2].find("start_codon") != -1 and words[6] == "+":
                genId += 1
            if words[2].find("stop_codon") != -1 and words[6] == "-":
                genId += 1
            if int(words[3]) >= startindex and int(words[3]) <= endindex:
                words[0] += "_%d" % genId
                words[1] = program
                words[3] = str(int(words[3]) - startindex)
                words[4] = str(int(words[4]) - startindex)
                print >>outfile, "\t".join(words)
            if int(words[3]) > endindex:
                break

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--startindex', help="Start index of the part to extract. Entry Indices will be adjusted to this value, meaning, here you should be precise. Take the value: Sbjct_Index - Query_Index", type=int, nargs="?")
    parser.add_argument('-l', '--length', help="Start  index of the part to extract. Entry Indices will  be adjusted to this  value, meaning, here you should be precise. Take  the value: Sbjct_Index -  Query_Index", type=int, nargs="?")
    parser.add_argument('-e', '--endindex', help="End index. Only Entries smaller than this value are included", type=int, nargs="?")
    parser.add_argument('-p', '--program', help='The name of the program which generated the GTF file, e.g. twinscan or CONTRAST',required=True)
    parser.add_argument('infile', help="The GTF input file.",)
    parser.add_argument('outfile', help="The GFF output file.", nargs="?")
    args = parser.parse_args()
    #Check argument consistency
    numLengthArgs = (1 if args.startindex else 0) + (1 if args.endindex else 0) + (1 if args.length else 0)
    if numLengthArgs < 2:
      parser.print_help()
      print "You need to specify at least two of --startindex, --length and --endindex"
      sys.exit(1)
    if args.startindex is not None and args.endindex is not None and args.startindex > args.endindex:
        parser.print_help()
        print('Check your start and end indices.')
        sys.exit(1)
    if args.length is not None and args.length < 1:
        parser.print_help()
        print('Length too short')
        sys.exit(1)
    if args.length is not None and args.startindex is not None and args.endindex is not None and (args.endindex - args.startindex) != args.length:
        parser.print_help()
        print('Length does not match start/end index.')
        sys.exit(1)
    if args.startindex is None: args.startindex = args.endindex - args.length
    if args.endindex is None: args.endindex = args.startindex + args.length
    #Execute the converter
    gtf2gff(args.infile, args.outfile, args.startindex, args.endindex, args.program)
