#!/usr/bin/env python3

import click
import os
import sys
from Bio import SeqIO
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command()
@click.version_option(version='1.0.0')
@click.option('-i', '--input', default='all', type=str, hide_input=False, nargs=1, required=True, help="Fasta input")
@click.option('-o', '--output', default='all', type=str, hide_input=False, nargs=1, required=False, help="Fasta Output")
@click.option('-w', '--width', default='all', type=int, hide_input=False, nargs=1, required=False, help="Fasta Output")

#type=click.Choice(['md5', 'sha1'])
#type=click.IntRange(0, 20, clamp=True)
#type=click.File('rb')
#type=click.File('wb')




def fastaFormatter (input, width, output):
	
	click.echo('input fasta: %s' % input)
	click.echo('input width: %i' % width)
	click.echo('Output fasta: %s' % output)
	outFa=open (output, "w+")
	with open(input, "r+") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			outFa.write (">" + record.description + os.linesep)
#			print (">" + record.description)
#			print (repr(record.seq))
#			type(record.seq)
			if width ==0:
				outFa.write(str(record.seq) + os.linesep)
			elif width >0 :
				startSeq=0
				while startSeq<len(record.seq):
					endSeq=startSeq+width
					if endSeq>len(record.seq):
						endSeq=len(record.seq)
					fastaStr=record.seq[startSeq:endSeq]
					
#					print (fastaStr)
					outFa.write(str(fastaStr) + os.linesep)
					startSeq+=width
			else:
				sys.stderr.write("Error: invalid width parameter: " + width)
	sys.stderr.write("Info: gracefully done" + os.linesep)



if __name__=='__main__':
	fastaFormatter()
