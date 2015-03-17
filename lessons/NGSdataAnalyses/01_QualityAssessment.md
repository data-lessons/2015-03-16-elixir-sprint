---
layout: page
title: Quality Assessment
subtitle: Checking quality sequence data
minutes: 30
---
> ## Learning Objectives {.objectives}
>
> * Understanding a FASTQ file structure
> * Controlling the Quality and encoding of a FASTQ file
> * Switching to one Quality encoding type to another one

#Picking up Data

Download the [dataset][dataLink] and uncompress it using
~~~{.bash}
~$ unzip data.zip
~~~
The original FASTQ data are in the *Data* folder (*all_seq_1.fastq* and *all_seq_2.fastq*).

#Understanding FASTQ Format

Take a look to the FASTQ files content using the *head* or *tail* command utilities

~~~{.raw}
~$ head data/all_seq_1.fastq
@RC10_HWUSI-EAS454_0006:1:99:16639:1487#TAGCTT/1
TTCTTGTGTAGATTGGGAAATTTCAGTTGGACTGCATCAATGGGGATCCCCTAGTTGGCCTCAGCAAGTGTGGAAG
+
fffdfefffafefffffffcfeffdfffffefefffefeefbefef`effeeefeff`fc^eeeea`d`dadbbad
@RC10_HWUSI-EAS454_0006:1:110:8483:17695#TAGCTT/1
GAAAAGTGCATCCACACTTGCTGAGGCCAACTAGGGGATCCCCATTGATGCAGTCCAACTGAAATTTCCCAATCTA
+
fffffefcffffffffffffffffffffffedeedfeefeffff^feefffdffcedfddfdfdeddcfea\dade
@RC10_HWUSI-EAS454_0006:1:25:5030:7161#TAGCTT/1
NTCAATTCTTGTGTAGATTGGGAAATTTCAGTTGGACTGCATCAATGGGGATCCCCTAGTTGGCCTCAGCAAGTGA
~~~

It will automatically display in the terminal the 10 first lines of the *all_seq_1.fastq* file. Using *tail* will display the last 10 lines instead.

Individual reads are represented on four lines:
> 1- Name of the read preceeded by a "@"
> 2- Sequence of the read
> 3- A line with '+', that can be followed by the name of the read again
> 4- The Quality score of each base encoded in ASCII format

The Quality is encoded in ASCII format, as followed (from https://en.wikipedia.org/wiki/FASTQ_format#Encoding)

![FASTQ Quality Encoding](../../img/NGSmapping_fastqEncoding.png  "FASTQ Quality Encoding")

The Quality value (from 0 to 40 generally) represents the probability that the given base is not an erroneous one. The exact formula is Q<sub>sanger</sub> = -10 log<sub>10</sub> *p*

#Quality control using FASTQC

FASTQC can be used to obtain different informations on the Fastq files: read number, size, quality per base and per sequence, GC content per base and per sequence, N content, overrepresented sequences and k-mer content. It can be used either through grapical way or in command line. We will use here the command line to obtain those informations.

~~~{.bash}
~$ fastqc -o 1_fastqc/ all_seq_1.fastq
~~~

![FASTQC_gui_output](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc.png  "FASTQC GUI Output")

Here the software will provide the different informations on the *all_seq_1.fastq* file and will report those informations in the *1_fastqc* folder (that will be created) in an HTML file.

> ## Challenge {.challenge}
>
> 1. What is the type of encoding for the quality ?
> 2. Which are the overrepresented sequences and their origin ?
> 3. How many sequences are present in the *all_seq_2.fastq* file ? You can also use the following bash command to obtain this information.
~~~{.bash}
~$ for i in `wc -l data/all_seq_1.fastq |cut -f1 -d" "` ;do echo $(($i/4));done
~~~

# Changing the Quality encoding type

As you can observe in the FASTQ files and in the FASTQC report, the files are encoded in PHRED+64. We need to change the encoding type before launching the cleaning step.

We will use the *seqret* tool from the *EMBOSS* suite to switch from PHRED+64 to PHRED+33.

~~~{.bash}
~$ seqret fastq-illumina::all_seq_1.fastq fastq-sanger::all_seq_1_sanger.fastq
~~~

> ## Challenges {.challenge}
>
> 1. Perform the same changes on the *all_seq_2.fastq* file
> 2. Repeat the FASTQC experiment and confirm the PHRED+33 encoding type

[dataLink]:../../data/biology/NGSdata/data.zip
