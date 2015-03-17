---
layout: page
title: SNP Calling
subtitle: SNP Calling using mpileup utility
minutes: 20
---
# What is SNP calling ?

A SNP (Single Nucleotide Polymorphism) is a single base change in an individual (or more) compared to a reference sequence.

SNP calling is the informatics process aimed at detecting such mutations and at outputting them in a human-readable format. The process uses a BAM file and a reference, and numerous tools are able to call SNPs, but with different methods behind them.

## Why *mpileup* ?

We choose the *SAMtools mpileup* utility because it is a quite fast tool in the frame of the lesson. However, more efficient and precise tools exist for SNP calling (*GATK*, *MrBayes*, ...).

##Running *mpileup*

The mpileup command is as follow:

~~~{.bash}
~$ mkdir 5_SNPcalling
~$ samtools mpileup -u -f data/reference.fas 4_SAMtools/rmdupProperlyMapped.bam > 5_SNPcalling/rawSNP.bcf
~~~
The output file is a bcf (for Binary Call Format), a compressed binary file. To be read by humans, it has to be transformed in VCF (for [Variant Call Format][vcfLink]) using the *bcftools view* utility:

~~~{.bash}
~$ bcftools view ­‐v ­‐c ­‐g 5_SNPcalling/rawSNP.bcf > 5_SNPcalling/rawSNP.vcf
~~~
A classical VCF file looks as follow:

![VCFpicture](http://bioinf.comav.upv.es/courses/sequence_analysis/_images/vcf_format.png "VCF Picture")

After the header lines, every line will represent a position with a SNP, compared to the reference.

#BE CAREFUL!!

The SNP calling will provide a list of **putative** SNPs, that may be true or not. Post-treatment are requested to go further on the biological relevance of the detected SNPs.

> ## Challenges {.challenge}
>
> 1. Count the number of SNPs using the *grep* and *wc* commands
> 2. Can you find heterozygous SNP ? Are they expected for a 99% autogamous plant ?

[vcfLink]: http://samtools.github.io/hts-specs/VCFv4.2.pdf
