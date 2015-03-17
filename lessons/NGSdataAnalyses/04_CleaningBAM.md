---
layout: page
title: CleaningBAM
subtitle: Obtaining a BAM ready for SNP calling
minutes: 10
---
> ## Learning Objectives {.objectives}
>
> * Understand that the already obtained BAM contains biases
> * To be able to clean a BAM file for SNP calling

# Removing the *not properly mapped* in pair reads
The mapping step has provide a SAM/BAM file which contains all the reads provided, whatever they are mapped or not, and properly mapped according to reference and options or not.
Here we will remove all read pairs that are not properly mapped **in pair**, meaning that:
> * the two mates of the given pair respect the mapping constraints,
> * that their relative orientation and distance of their outermost bases is correct.

For that we will re-use the *SAMtools view* utility but with different options.

~~~{.bash}
~$ samtools view -o 4_SAMtools/properlyMapped.bam -f 0x02 4_SAMtools/sorted.bam
~~~

> ## About the "-f 0x02" option {.callout}
>
> The *-f 0x02* option means that only reads with the *0x02* flag will be retained. This flag (see [SAM format][samSpecLink] for more description) means that the read and its mate are **properly mapped in pair**.

# Removing duplicated reads

Reads (or pairs) are considered as duplicates if their coordinates are exactly the same once mapped on the reference sequence. Such events arise because from different origins, that will not be discussed more here (see ####REF### for more informations).

We will use the *SAMtools rmdup* utility for that

~~~{.bash}
~$ samtools rmdup 4_SAMtools/properlyMapped.bam 4_SAMtools/rmdupProperlyMapped.bam
~~~

The file is now ready to be used for SNP calling.

[samSpecLink]: http://samtools.github.io/hts-specs/SAMv1.pdf
