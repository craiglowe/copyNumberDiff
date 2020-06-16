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

This code uses the Kent libraries from UCSC as well as the GNU Scientific Library and R, so those must be installed
to compile copyNumberDiff and its helper programs.

<ol>
<li> Download and compile the Kent Libraries

<ol>
<li> Set environment variables

<ol>
<li> Check the setting of your machtype variable with:<br />
echo $MACHTYPE<br />
it should be something in this list:<br />
i386 i686 sparc alpha x86_64 ppc<br />
If it is not, set your machtype variable.
If you are not sure what your MACHTYPE variable should be, you can try typing "uname -m", "uname -p", or "uname -a"
for some hints and pick the element from the previous list that most closely matches the output.  For most people
it will be "x86_64".
</ol>

<li> Go to a folder on your computer where you want the kent source tree to reside and type:<br />
git clone git://genome-source.cse.ucsc.edu/kent.git<br />
to download the repository onto your own computer.
<li> go to the src directory within the kent source repo that you just cloned:<br />
cd kent/src<br />
<li> Compile the libraries<br />
make libs
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
<li> Unpack the file you downloaded and go into the directory
<li> Configure the make files and you may want to use the prefix option to change the install path:<br />
./configure --prefix=/home/lowec/src/gsl/gsl-1.16_install
<li> Compile the source:<br />
make
<li> Move the library to the install location:<br />
make install
</ol>

<li> Download and compile R

<ol>
<li> Get R:<br />
wget http://cran.r-project.org/src/base/R-3/R-3.2.3.tar.gz
<li> Unpack the file you downloaded and go into the directory
<li> ./configure
<li> cd src/nmath/standalone
<li> make static
</ol>

<li> Compile copyNumberDiff and helper programs
<ol>
<li> Edit the copyNumberDiff makefile so that it points to the kent and gsl libraries on your system.  These
are the six lines you will have to modify:<br />
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a
HG_INC += -I/home/lowec/src/gsl/gsl-1.16_install/include
L += /home/lowec/src/gsl/gsl-1.16_install/lib/libgsl.a /home/lowec/src/gsl/gsl-1.16_install/lib/libgslcblas.a
HG_INC += -I/home/lowec/src/R/R-3.1.2/include
L += /home/lowec/src/R/R-3.1.2/src/nmath/standalone/libRmath.a

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

<li> bamToGcStats is not needed for the workflow, but it can show you a summary
of the GC content reads from a bam file.  This can give you the GC-bias of the library
when compared against the stats you just calculated of the reference genome.
<ol>
<li> in.bam - the bam file to analyze.
<li> readLength - sequencing read length of library in bam file.
<li> out.stats - this is a text output file summarizing the sequencing library
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

<li> wgHmm is in a separate repository since it is a general purpose HMM/transducer
program for genome-wide analyses.  Given the emission probabilities from above
and a set of transition probabilities (penalties), it will report the most likely
series of states (Viterbi path).  This will be the final output that will report
the regions of repeatedly different copy number between the two groups.
<ol>
<li> -inputIsLog - the transition matrix and the emisison probabilities are
interpreted as being in log-space.  This is probably a good idea in most cases.
<li> chrom.sizes - tab separated file where the first column is the name
of the chromosome and the second column is the length
<li> noGaps.bed - bed file of regions no overlapping an assembly gap
<li> numberOfStates - number of states (25 in this case)
<li> transition.matrix - is a tab separated file where column X row Y gives the
transition probabiliy from state X to state Y
<ol>
<li> setting these parameters is the hardest part of running this program.  I do
not know of a good way to quickly calculate them.  In the past I have set them by
simulating data sets where I know the location of repeated copy number variation.
I then make the transition probabilities low enough to reduce false positives, while
maintaining the ability to detect the true positives.  These parameters may go as low
as -600 (log-space).
</ol>
<li> emissionProbs.wig - based on a fixed-step wiggle file, but with multiple columns
one for each state.  This is the output of copyNumberDiff
<li> output.bed - these are the final calls of the program and pipeline.  The first
three columns give the genome coordinates and the next column is a number that
describes the copy number.  The number is equal to (5 * copiesGroup2 + copiesGroup1).
The copies are 0 for homozygous deletion, 1 for heterozygous deletion, 2 for reference,
3 for heterozygous duplication, and 4 for homozygous duplication.  For example:
14 = 5 * 2 + 4 = reference for group2, duplication for group1.
</ol>

</ol>


References
==========

This code has not yet been published.

