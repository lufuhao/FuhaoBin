#!/usr/bin/env perl
use strict;
use warnings;
use constant LICKHELP =><<EOM;
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

#装饰一个函数，使之成为命令行借口
#type=click.Choice(['md5', 'sha1'])
#type=click.IntRange(0, 20, clamp=True)
#type=click.File('rb')
#type=click.File('wb')
\@click.command()
#添加命令行选项
\@click.version_option(version='1.0.0')
\@click.option('--file', default='all', type=string/int/float, hide_input=False, prompt="Your choice: ", nargs=1, required=True, help="input")

def hello (count, name):
	"""Simple program that greets NAME for a total of COUNT times."""
	for x in range(count):
		click.echo('Hello \%s!' % name)



if __name__=='__main__':
	hello()
EOM
print LICKHELP;
