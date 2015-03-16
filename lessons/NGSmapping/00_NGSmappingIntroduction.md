---
layout: page
title: NGS Data Analysis 
subtitle: Main presentation
minutes: 10
---
> ## Learning Objectives {.objectives}
>
> * Understanding a FASTQ file structure
> * Quality assessement and cleaning data
> * Mapping NGS data upon a reference genome

# NGS Mapping

Typically, NGS data treatments are structured in 3 main steps:
1. Quality Assessment
2. Mapping upon a Reference
3. FInal analyses

The two first steps are standard operations, whereas the last one is specific to scientific questions.

The goal of this practical is to learn how to (*i*) perform basic Quality assessment using FASTQC and CutAdapt, then (*ii*) to be able to perform a mapping using bwa upon a reference sequence, and (*iii*) then just extracting the SNP using SAMtools as an example of third step.

#Before we get started
You need to have access to the following tools:

> * Basic Shell utilities: unzip, tar, cp, mkdir, mv, wc, grep, echo, rm, head, tail, less, cat, sed
> * FASTQC
>> For Linux or windows (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip)
>> For Mac (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.dmg)
> * CutAdapt (https://github.com/marcelm/cutadapt/archive/master.zip)
> * bwa (http://sourceforge.net/projects/bio-bwa/files/latest/download?source=files)
> * Re-pairing tool (http://***/pairing.pl)
> * EMBOSS packages (ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz)
> * SAMtools (http://sourceforge.net/projects/samtools/files/latest/download?source=files)

You must have also a suitable Java JRE 1.7xx or more.

[FASTQC] [fastqcLink] is a quality control tool for high throughput sequence data, working in either a graphical way or in command line. It will provide a HTML output giving informations about quality encoding, data quality, overrepresented sequences and so on.

[CutAdapt][cutadaptLink] is a tool aimed at removing adaptors sequences and at quality-based trimming of reads. 

[bwa][bwaLink] is a tool for mapping short as long reads upon a reference sequence using a Burrows-Wheeler Transform method.

[pairing.pl][pairingLink] is a simple Perl script able to re-pairing NGS sequences after the cleaning step.

[EMBOSS][embossLink] is a suite of basic tools for sequence manipulations.

[SAMtools][samLink] is a suite of small program able to work with SAM and BAM files (mapping files): converting SAM to BAM, extraction of a part of SAM/BAM, SNP extraction through pileup...


#More resources
> * Links to articles
> * Websites

[fastqcLink]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[cutadaptLink]: https://code.google.com/p/cutadapt/
[bwaLink]: http://bio-bwa.sourceforge.net/
[embossLink]: http://emboss.sourceforge.net/
[pairingLink]: http://***/pairing.pl
[samLink]: http://www.htslib.org/
