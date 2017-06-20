#!/usr/bin/env perl

use Getopt::Long;
use warnings FATAL => 'all';
use strict;

#options
my $optVerbose;

sub usage
{
    print STDERR "
usage:   printTransMatrix.pl deletionPenalty duplicationPenalty expectedLength outputFile
options:
     -verbose
";
    exit(1);
}

sub checkOptions
{
    # Make sure command line options are valid/supported.
    my $ok = GetOptions("verbose" => \$optVerbose);
    if(! $ok){&usage();}
    if(!(defined $optVerbose)){$optVerbose = 0;}

    if($optVerbose){print(STDERR "Options look ok\n");}
}

sub abort
{
    my ($message) = @_;
    print(STDERR "$message\n");
    exit(1);
}

sub makeMatrix
{
	my ($delPen, $dupPen, $expectedLength, $outputFile) = @_;

	my ($i, $j, $x, $y, $penalty, $selfTrans, $backTrans);

	$selfTrans = log(1.0 - 1.0/$expectedLength);
	$backTrans = log(1.0/$expectedLength);
	my @pens = ($delPen, 0, $dupPen);

	open(Fout,">$outputFile");
	
	for($x=0; $x <= $#pens; $x++)
	{
		foreach($y=0; $y <= $#pens; $y++)
		{
			if($x != $y){$penalty = $pens[$x] + $pens[$y];}
			else{$penalty = $pens[$x];}
			if($x == $#pens && $y == $#pens){print(Fout "$penalty\n");}
			else{print(Fout "$penalty\t");}
		}
	}
	
	for($i=0; $i <= $#pens; $i++)
	{
		for($j=0; $j <= $#pens; $j++)
		{
			for($x=0; $x <= $#pens; $x++)
			{
				for($y=0; $y <= $#pens; $y++)
				{
					if($i == 1 && $j == 1 && $x == 1 && $y == 1){print(Fout "0.0");}
					elsif($i == $x && $j == $y){print(Fout "$selfTrans");}
					elsif($x == 1 && $y == 1){print(Fout "$backTrans");}
					else
					{
						if($x != $y){$penalty = $pens[$x] + $pens[$y];}
						else{$penalty = $pens[$x];}
						print(Fout "$penalty");
					}

					if($x == $#pens && $y == $#pens){print(Fout "\n");}
					else{print(Fout "\t");}
				}
			}
		}
	}
	close(Fout);
}

######################################################
# main
######################################################


&checkOptions();

if (($#ARGV + 1) != 4){&usage();}
my $delPenalty = $ARGV[0];
my $dupPenalty = $ARGV[1];
my $expectedLen = $ARGV[2];
my $outputFilename = $ARGV[3];

&makeMatrix(-1*$delPenalty, -1*$dupPenalty, $expectedLen, $outputFilename);

