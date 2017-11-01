#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use DBI;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20141016

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--db    <database_name>
		Default
	--verbose
		Detailed output for trouble-shooting;
	--version|v!
		Print current SCRIPT version;

Example:
	perl $0 

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $dbname, $host, $port, $username, $password);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"db:s" => \$dbname,
	"host:s" => \$host,
	"port:i" => \$port,
	"user:s" => \$username,
	"password:s" => \$password,
#	!:s:i
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
$dbname='b2gdb' unless(defined $dbname);
$host="localhost" unless (defined $host);
$port=3306 unless (defined $port);
$username='blast2go' unless (defined $username);
$password='blast4it' unless (defined $password);



### input and output ################################################



### Main ############################################################
# Connect to the database.
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host:$port", $username, $password, {'RaiseError' => 1});
$dbh->do("drop tabletruncate table if exists `gi2uniprot`");

 if exists `gene2accession`;


  # Drop table 'foo'. This may fail, if 'foo' doesn't exist
  # Thus we put an eval around it.
  eval { $dbh->do("DROP TABLE foo") };
  print "Dropping foo failed: $@\n" if $@;

  # Create a new table 'foo'. This must not fail, thus we don't
  # catch errors.
  $dbh->do("CREATE TABLE foo (id INTEGER, name VARCHAR(20))");

  # INSERT some data into 'foo'. We are using $dbh->quote() for
  # quoting the name.
  $dbh->do("INSERT INTO foo VALUES (1, " . $dbh->quote("Tim") . ")");

  # same thing, but using placeholders (recommended!)
  $dbh->do("INSERT INTO foo VALUES (?, ?)", undef, 2, "Jochen");

  # now retrieve data from the table.
  my $sth = $dbh->prepare("SELECT * FROM foo");
  $sth->execute();
  while (my $ref = $sth->fetchrow_hashref()) {
    print "Found a row: id = $ref->{'id'}, name = $ref->{'name'}\n";
  }
  $sth->finish();

  # Disconnect from the database.
  $dbh->disconnect();



#####################################################################
###                         sub functions                         ###
#####################################################################
### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependency:
###Note:



#idmapping.tb
#This table includes the following IDs (or ACs) delimited by tab:
#1. UniProtKB accession               A0A1L2CPP7
#2. UniProtKB ID                      A0A1L2CPP7_9BIVA
#3. EntrezGene	
#4. RefSeq
#5. NCBI GI number                    1000001402
#6. PDB
#7. Pfam
#8. GO
#9. PIRSF
#10. IPI
#11. UniRef100
#12. UniRef90
#13. UniRef50
#14. UniParc
#15. PIR-PSD accession
#16. NCBI taxonomy
#17. MIM
#18. UniGene
#19. Ensembl
#20. PubMed ID
#21. EMBL/GenBank/DDBJ
#22. EMBL protein_id

