#!/bin/csh -ef

# edit the below lines to point to the programs on your system
set evolveGenome = "/home/cl454/src/useKent/evolveGenome/evolveGenome"
set countKmers = "/home/cl454/src/useKent/copyNumberDiff/countKmers"
set faToGcStats = "/home/cl454/src/useKent/copyNumberDiff/faToGcStats"
set copyNumberDiff = "/home/cl454/src/useKent/copyNumberDiff/copyNumberDiff"
set wgHmm = "/home/cl454/src/useKent/wgHmm/wgHmm"
set art = "/home/cl454/src/art_src_MountRainier_Linux/art_illumina"
set bwa = "/home/cl454/src/bwa-0.7.17/bwa"
set samtools = "/home/cl454/src/samtools-1.10/samtools"
set overlapSelect = "/home/cl454/bin/x86_64/overlapSelect"

# these files should already be in the folder.  They are one
# chromosome from the stickleback assembly
set chromSizes = "gasAcu1.chrXXI.chromSizes"
set startingGenome = "gasAcu1.chrXXI.fa"
set noGap = "gasAcu1.chrXXI.noGap.bed"

# this is how you control how many individuals are in each group
set groupOneRef = "ref/g1.1.bam,ref/g1.2.bam,ref/g1.3.bam,ref/g1.4.bam,ref/g1.5.bam"
set groupTwoRef = "ref/g2.1.bam,ref/g2.2.bam,ref/g2.3.bam,ref/g2.4.bam,ref/g2.5.bam"
set groupTwoDel = "del/g2.1.bam,del/g2.2.bam,del/g2.3.bam,del/g2.4.bam,del/g2.5.bam"
set groupTwoDup = "dup/g2.1.bam,dup/g2.2.bam,dup/g2.3.bam,dup/g2.4.bam,dup/g2.5.bam"

echo "...Placing Deletions and Duplications In Genomes..."
$evolveGenome $startingGenome $noGap genomeWithDeletions.fa deletionLocations.bed -dels=100 -minDel=30 -maxDel=30
$evolveGenome $startingGenome $noGap genomeWithDuplications.fa duplicationLocations.bed -dups=100 -minDup=200 -maxDup=200

echo "...Counting Kmers..."
$countKmers -verbose=2 36 $startingGenome $noGap stdout | gzip > 36mers.gz

echo "...Analyzing GC content..."
$faToGcStats -verbose=2 $startingGenome $noGap 36 gc.log -wiggle=/dev/stdout | gzip > 36gc.gz

echo "...Building BWA index..."
$bwa index $startingGenome

echo "...Simulating Reads..."
mkdir -p ref dup del
foreach f (`echo $groupOneRef | tr ',' ' '`)
	set n = $f:t:r
	echo $n
	$art -i $startingGenome -l 36 -f 3 -o ref/$n -na
	$bwa mem $startingGenome ref/$n.fq | $samtools view -F 4 -bS - > $f
end
foreach f (`echo $groupTwoRef | tr ',' ' '`)
	set n = $f:t:r
	echo $n
	$art -i $startingGenome -l 36 -f 3 -o ref/$n -na
	$bwa mem $startingGenome ref/$n.fq | $samtools view -F 4 -bS - > $f
	$art -i genomeWithDeletions.fa -l 36 -f 3 -o del/$n -na
	$bwa mem $startingGenome del/$n.fq | $samtools view -F 4 -bS - > del/$n.bam
	$art -i genomeWithDuplications.fa -l 36 -f 3 -o dup/$n -na
	$bwa mem $startingGenome dup/$n.fq | $samtools view -F 4 -bS - > dup/$n.bam
end

echo "...Calculating Emission Probabilities..."
$copyNumberDiff $noGap 36mers.gz 36gc.gz $groupOneRef $groupTwoRef stdout -noHet -verbose=3 -gcPseudo=10 -logProbs | gzip > probs.ref.noHet.wig.gz
$copyNumberDiff $noGap 36mers.gz 36gc.gz $groupOneRef $groupTwoDel stdout -noHet -verbose=3 -gcPseudo=10 -logProbs | gzip > probs.del.noHet.wig.gz
$copyNumberDiff $noGap 36mers.gz 36gc.gz $groupOneRef $groupTwoDup stdout -noHet -verbose=3 -gcPseudo=10 -logProbs | gzip > probs.dup.noHet.wig.gz

echo "...Running Transducer..."
set dels = 1
set dups = 1
set delPen = 200
set dupPen = 200
while ( $dels > 0 && $dups > 0 )
	echo "Deletion and duplication penalties are: $delPen and $dupPen"
	./printTransMatrix.pl $delPen $dupPen 50 current.mat
	$wgHmm -inputIsLog $chromSizes $noGap 9 current.mat probs.ref.noHet.wig.gz output.bed
	cat output.bed \
	| perl -ne 'chomp($_);@w=split("\t",$_);if(($w[3] eq "1")||($w[3] eq "2")||($w[3] eq "3")||($w[3] eq "6")){print("$_\n");}' > output.dels.bed
	cat output.bed \
	| perl -ne 'chomp($_);@w=split("\t",$_);if(($w[3] eq "2")||($w[3] eq "5")||($w[3] eq "6")||($w[3] eq "7")){print("$_\n");}' > output.dups.bed
	set dels = `cat output.dels.bed | wc -l`
	set dups = `cat output.dups.bed | wc -l`
	echo "There were $dels deletions and $dups duplications"
	if ( $dels > $dups ) then
		@ delPen = $delPen + 50
	else if ( $dups > $dels ) then
		@ dupPen = $dupPen + 50
	else if ( $dups > 0 ) then
		@ dupPen = $dupPen + 50
	endif
end

echo "The final deletion and duplication penalties are: $delPen and $dupPen"

./printTransMatrix.pl $delPen $dupPen 50 $delPen.$dupPen.mat
foreach mut (del dup)
	$wgHmm -inputIsLog $chromSizes $noGap 9 $delPen.$dupPen.mat probs.$mut.noHet.wig.gz output.$mut.bed
	cat output.$mut.bed | perl -ne 'chomp($_);@w=split("\t",$_);if(($w[3] ne "0")&&($w[3] ne "4")&&($w[3] ne "8")){print("$_\n");}' > output.$mut.calls.bed
end

set delNumer = `$overlapSelect output.del.calls.bed deletionLocations.bed stdout | wc -l`
set delDenom = `cat deletionLocations.bed | wc -l`
set dupNumer = `$overlapSelect output.dup.calls.bed duplicationLocations.bed stdout | wc -l`
set dupDenom = `cat duplicationLocations.bed | wc -l`

echo "The program found $delNumer of $delDenom deletions and $dupNumer of $dupDenom duplications"

echo DONE

