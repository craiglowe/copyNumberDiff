copyNumberDiff
================

A genomics tool that lets you compare two groups of individuals to see if there are copy number differences
the correlate with group assignment.  These groups can be defined based on an ecological variable, such
as fish that live in marine water, and fish that live in fresh water.  The groups could also be defined
based on the presence and absense of a disease.  In the case of tumor and matched normal samples, the
groups would not even need to be made up of individuals, but could be made up of different tissue samples
from one or more individuals.


Code I Did Not Write
============
One of the helper programs uses a Bloom filter.  I used the MurmurHash3 family of hash functions written
by Austin Appleby while implementing the bloom filter and have included his murmur3.h and murmur3.c
files in this repository.  Austin has placed this work in the public domain.


Installation
============

This code uses the Kent libraries from UCSC as well as the GNU Scientific Library, so those must be installed
to compile copyNumberDiff and its helper programs.

<ol>
<li> Download and compile the Kent Libraries

<ol>
<li> Check the setting of your machtype variable with:<br />
echo $MACHTYPE<br />
it should be something in this list: i386 i686 sparc alpha x86_64 ppc.  If it is not, set your machtype variable.
<li> Go to a folder on your computer where you want the kent source tree to reside and type:<br />
git clone git://genome-source.cse.ucsc.edu/kent.git<br />
to download the repository onto your own computer.
<li> go to the src/lib directory within the kent source repo that you just cloned:<br />
cd kent/src/lib<br />
<li> Compile the libraries<br />
make
<li> If this was successful, you should have a file here:<br />
kent/src/lib/x86_64/jkweb.a<br />
the x86_64 will be the machtype of your machine.</br />
</ol>

If this was not successful then you should look at the build instructions in the kent repo itself
by looking at this file:<br />
kent/src/product/README.building.source

<li> Download and compile the GNU Scientific Library

<ol>
<li> Get GSL:<br />
wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
<li> Configure the make files and you may want to use the prefix option to change the install path:<br />
./configure --prefix=/home/lowec/src/gsl/gsl-1.16_install
<li> Compile the source:<br />
make
<li> Move the library to the install location:<br />
make install
</ol>

<li> Compile copyNumberDiff and helper programs
<ol>
<li> Edit the copyNumberDiff makefile so that it points to the kent and gsl libraries on your system.  These
are the four lines you will have to modify:<br />
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc<br />
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a<br />
HG_INC += -I/home/lowec/src/gsl/gsl-1.16_install/include<br />
L += /home/lowec/src/gsl/gsl-1.16_install/lib/libgsl.a /home/lowec/src/gsl/gsl-1.16_install/lib/libgslcblas.a

<li> Compile copyNumberDiff and helper programs:<br />
make all

<li> Running copyNumberDiff or any of the helper programs with no parameters should give a brief help message.
</ol>
</ol>


Usage overview
==========
<ol>
<li> Starting point
<ol>
<li> Reference genome
<li> Files of reads mapped to the reference genome
<li> Read length used during sequencing
<ol>
<li> This is the length of a read.  If you are using Illumina technology, then it will probably be: 36,
76, 100, 150, etc.  At this point the program does not take into account paired-end or mate-pair technology
and how they aids in mapping certainty.  For 100bp single-end reads, or 100bp paired-end reads, the
read length is 100.  At this point the program assumes that all libraries where sequenced with the same
length.  Accomodating mixed read lengths would be possible, but this has not been implemented yet.
</ol>
</ol>

<li> We will use the countKmers program to help estimate the mappablity of
each base in the reference genome.  Running the program with no arguments
will display a help message.
<ol>
<li> kmerLength - this is the read length
<li> in.fa - this is the reference genome in fasta format
<li> noGap.bed - this is a bed file of all regions in the
genome that are not assembly gaps.
<li> kmers.wig - this is the output file.  You can use
stdout as the output filename to send it to standard out.
Doing this and piping it through gzip might be a good idea
because the file will usually be large.
</ol>

<li> We will use the faToGcStats program to estimte how likely the coverage
of a base is to be affected by GC-bias in the sequencing library.
<ol>
<li> in.fa - this is the reference genome in fasta format
<li> noGap.bed - this is a bed file of all regions in the
genome that are not assembly gaps.
<li> windowLength - this is the read length
<li> out.stats - this is a text output file summarizing the whole genome
<li> -wiggle (filename) - this is the output file that will go on to be
used by copyNumberDiff.
</ol>

<li> copyNumberDiff will use the results of the two previous program, along with
the bam files of each individual, to compute the emission probablilities of the transducer.
<ol>
<li> noGap.bed - this is a bed file of all regions in the
genome that are not assembly gaps.
<li> kmerMismaps.wig - this is the output from countKmers.  It can be compressed with
gzip and the program will figure this out if the filename ends with .gz.
<li> gcContent.wig - this is the output from faToGcStats.  It can be comprssed with
gzip and the program will figure this out if the filename ends with .gz.
<li> group1_1.bam,...,group1_N.bam - this is a comma separated list of all the bam
files in group 1.
<li> group2_1.bam,...,group2_N.bam - this is a comma separated list of all the bam
files in group 2.
<li> outputEmissionProbs.wig - this is the output of the program and it desribes
the emission probabilities of the transducer.  It may be a good idea to send this
to stdout and pipe it through gzip.
<li> Good to use the -logProbs option to prevent repeated multiplications of probabilities
from going to zero.
</ol>

</ol>


References
==========

This code has not yet been published.

