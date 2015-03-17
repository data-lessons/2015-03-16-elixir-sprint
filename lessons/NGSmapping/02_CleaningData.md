---
layout: page
title: Cleaning Data
subtitle: Cleaning sequence data using CutAdapt
minutes: 10
---
> ## Learning Objectives {.objectives}
>
> * Understanding the cleaning methods and goals
> * Learning to use CutAdapt

#Why Cleaning ?

Raw sequence data are still contaminated by various bias: presence of remaining sequencing adaptors, low quality bases and so on. To avoid an impact in the final analyses, we need to remove those artifacts and to control the quality of the remaining reads that will be analyzed forward.

#How to clean ?

A lot of tools are freely available for cleaning NGS data (FastX-toolkit, *Trimmomatic*, NGS QC toolkit, etc...). Here, we will use CutAdapt, a tool originally aimed to remove adaptor sequences but that also can trim the low-quality bases.

##The way CutAdapt works

[CutAdapt][cutadaptLink] will first perform a scan of the sequence in order to identify homologous sequences to the list of adaptors provided by the user. The sequence of the adaptors will be removed then. CutAdapt will also trim the scanned sequence in order to retain ONLY bases with a Quality score higher than a given threshold.

## Use of CutAdapt

Here we will use the Solexa_mRNA_primers.txt file to generate the list of adaptors to provide to CutAdapt.

1. First, create the "2_Cutadapt" folder using the *mkdir* command
~~~{.bash}
~$ mkdir 2_Cutadapt
~~~
2. Generate the list of reverse complementary adaptors using the *revseq* tool from the *EMBOSS* suite
~~~{.bash}
~$ revseq -sequence Solexa_mRNA_primers.txt -outseq Solexa_mRNA_primers.txt.revers.compl -reverse Y -complement Y
~~~
3. Generate the **cutadapt.conf** configuration file using the *grep* and *sed* command utilities
~~~{.bash}
~$ grep -v "Solexa*" | sed -e 's/\n/\n-b/g' > cutadapt.conf
~~~
4. Add specific options to your current experiment to the **cutadapt.conf** file
~~~{.bash}
~$ echo -e "-O 7\n-m 20\n-q 20\n-e 0.1"
~~~

Now, the **cutadapt.conf** file is ready.

> ## Options explained {.callout}
>
> * -b: Finding this adaptor in 5' as in 3'
> * -O: minimum overlap between the adaptor and the read sequence
> * -m: Minimum size for a read sequence to be retained after trimming
> * -q: Quality score threshold
> * -e: Error rate for the adaptor sequence

## Running CutAdapt

Now the configuration file is ready, CutAdapt can be launched on the sequence file, using the following command for the first one:

~~~{.bash}
cutadapt $(<cutadapt.conf) all_seq_1_sanger.fastq > all_seq_1_sanger_trimmed.fastq
~~~

> ## Challenges {.challenge}
>
> 1. Repeat the CutAdapt step using the second file
> 2. You can check the modification of sequence data using FASTQC as in the previous [lesson 01_QualityAssessement][lessonO1NGS]

[cutadaptLink]: https://code.google.com/p/cutadapt/
[lessonO1NGS]: 01_QualityAssessment.md
