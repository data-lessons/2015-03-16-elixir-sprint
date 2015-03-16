---
layout: page
title: Quality Assessment
subtitle: Checking and cleaning data
minutes: 30
---
> ## Learning Objectives {.objectives}
>
> * Understanding a FASTQ file structure
> * Controlling the Quality and encoding of a FASTQ file
> * Removing adapters and low-quality reads and bases

#Picking up Data

Download the [dataset][dataLink] and uncompress it using
~~~{.bash}
~$ unzip data.zip
~~~
The original FASTQ data are in the *Data* folder (*all_seq_1.fastq *and *all_seq_2.fastq*).

#Understanding FASTQ Format

Take a look to the FASTQ files content using the *head* or *tail* command utilities

~~~{.bash}
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

<FONT FACE="courier">
<span style="color: purple">SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS</span>.....................................................
   ..........................<span style="color: green">XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX</span>......................
   ...............................<span style="color: blue">IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII</span>......................
   .................................<span style="color: orange">JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ</span>......................
   <span style="color: brown">LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL</span>....................................................
   !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 |.........................|..|........|..............................|......................|
  33.......................59.64.......73............................104....................126
 <span style="color: purple">  0........................26.31.......40                                </span>
 <span style="color: green" >  .......................-5....0........9.............................40 </span>
 <span style="color: blue"  >    .............................0........9.............................40 </span>
 <span style="color: orange">  ................................3.....9.............................40 </span>
 <span style="color: brown" >  0.2......................26...31........41                              </span>
 
  <span style="color: purple">S - Sanger        Phred+33,  raw reads typically (0, 40)</span>
  <span style="color: green" >X - Solexa        Solexa+64, raw reads typically (-5, 40)</span>
  <span style="color: blue"  >I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)</span>
  <span style="color: orange">J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
      with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
      (Note: See discussion above).</span>
  <span style="color: brown" >L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)</span>
</FONT>

The Quality value (from 0 to 40 generally) represents the probability that the given base is not an erroneous one. The exact formula is $Q_\text{sanger} = -10 \, \log_{10} p$.

#Quality control using FASTQC

FASTQC can be used to obtain different informations on the Fastq files: read number, size, quality per base and per sequence, GC content per base and per sequence, N content, overrepresented sequences and k-mer content. It can be used either through grapical way or in command line. We will use here the command line to obtain those informations.

~~~{.bash}
~$ fastqc -o 1_fastqc/ all_seq_1.fastq
~~~

![FASTQC_gui_output](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc.png  "FASTQC GUI Output")

Here the software will provide the different informations on the *all_seq_1.fastq* file and will report those informations in the *1_fastqc* folder (that will be created) in an HTML file.




~~~ {.python}
some code:
    to be displayed
~~~
~~~ {.output}
output
from
program
~~~
~~~ {.error}
error reports from program (if any)
~~~

and possibly including:

> ## Callout Box {.callout}
>
> An aside of some kind.

> ## Challenge Title {.challenge}
>
> Description of a single challenge.
> There may be several challenges
> that make reference to [Challenge Title](01-one.html#challenge-title).


[dataLink]:http://***/data.zip