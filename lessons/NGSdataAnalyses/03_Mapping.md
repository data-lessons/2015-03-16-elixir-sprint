---
layout: page
title: Mapping
subtitle: Mapping FASTQ files upon a reference FASTA sequence
minutes: 20
---
> ## Learning Objectives {.objectives}
>
> * Understanding what is a mapping for NGS data
> * To be able to launch a mapping for NGS data using *bwa*

# What is a mapping for NGS data ?

Mapping is the process used to align the NGS reads upon a reference sequence. It is aimed to identify the different possible positions (if not unique) of any given NGS sequence on the reference. The mapping tools are numerous (*bwa*, *SOAP*, *bowtie*, *MAQ*, and many others) but all of them will provide a SAM file as output (Sequence/Alignment Map).

So, to run a mapping you need at least a reference sequence and a NGS sequence file. Here, we will use pair-ended data to perform the mapping.

# Running a mapping using *bwa aln*

*bwa* software will take into account mismatches and indels (insertions and deletions) for mapping, as well as the quality of the reads.

First, the reference has to be indexed to be used in the mapping process. Run the following command to index it:

~~~{.bash}
~$ bwa index data/reference.fas
~~~

The mapping operation itself will be split in two: independent mapping for forward and reverse files, then re-compiling the possible position for the pairs in a single SAM file.

The independent mapping can be performed using the following commands:

~~~{.bash}
~$ mkdir 3_bwa
~$ bwa aln -f 3_bwa/forward.sai data/reference.fas 2_Cutadapt/forward.fastq
~~~
> ## About the aln options {.callout}
>
> Here we did not change the mapping constraints for the *aln* part. However, there are many parameters that can be modified to optimize the mapping or to take into account the genetic distance between the reference and the individual that has been sequenced for instance. The most common are:
> * **Edit distance** (-n): the number of mismatches that are tolerated between the read and the reference.
> * **Gap opening** (-o): the number of gaps allowed between the read and the reference.
> * **Number of occurences of best hits** (-R): maximum number of alignments reported for a single read.
> * **Read Quality threshold for mapping** (-q): lower quality reads will be trimmed to 35 bases.

> ##Challenge {.challenge}
>
> Perform the mapping for the reverse file in the same way.

Once the independent mapping for the forward and reverse files are obtained, we can launch the mapping of the pairs themselves and extract the SAM file:

~~~{.bash}
~$ bwa sampe -f 3_bwa/all_seq.sam data/reference.fas 3_bwa/forward.sai 3_bwa/reverse.sai 2_Cutadapt/forward.fastq 2_Cutadapt/reverse.fastq
~~~

> ## About the sampe options {.callout}
>
> Here again we did not change the mapping constraints for the *sampe* part. However, there are also many parameters that can be changed. The most common are:
> * **Insert size** (-a): the size of the library, i.e. the distance between the outermost bases between the two mates of a pair. It depends of the sequencing process itself.
> * **Maximum hits to output for paired reads** (-n): only *n* positions/alignments will be reported in the final SAM file, if multiple positions exist for the considered pair.
> * **Maximum hits to output for discordant pairs** (-N): only *N* positions/alignments per mate will be reported in the final SAM file, for discordant mapping only (too large insert size, each mate on a different chromosome, abnormal position such as FF or RF or RR).

The resulting file *3_bwa/all_seq.sam* is in [SAM format][samSpecLink], a tabular text file in which the positions and variations for each reads compared to the reference are reported. Using *bwa sampe* there is only one line per read. Some tools may report more than a line per read, if many positions for this read exist.

# Converting SAM to BAM and re-ordering

The SAM file is a human-readable text file, but can be really huge (up to hundreds of Gb). To save space and to faster access to its content, we will transform it in a BAM format (for Binary/Alignment Map), similar to a compressed file (we will save more than 50% of drive space).

## Transforming a SAM in a BAM

We will use the *SAMtools view* utility to perform this step, and creating a BAM file in the *4_SAMtools* directory:

~~~{.bash}
~$ mkdir 4_SAMtools
~$ samtools view -bS -o 4_SAMtools/unsorted.bam 3_bwa/all_seq.sam
~~~

## Re-ordering the BAM file

We will use the *SAMtools sort* utility for that:

~~~{.bash}
~$ samtools sort -f 4_SAMtools/sorted.bam 4_SAMtools/unsorted.bam
~~~

The resulting BAM file is sorted according to the reference, i.e. starting from the read mapped the most closely possible to the 5' end of the first sequence in the reference to the read mapped the most closely to the 3' end of the last sequence in the reference.


[samSpecLink]: http://samtools.github.io/hts-specs/SAMv1.pdf
