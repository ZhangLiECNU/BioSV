## BioSV (Breakpoint-based identification of Structural Variation)

### Table of Contents
* [Overview](#overv)
* [Requirements and installation](#req)
* [Usage](#use)
* [Example](#exm)
* [Test data](#tes)


### <a name="overv"></a>Overview

Breakpoint-based identification of Structural Variation (BioSV), is an accurate and efficient SV caller, which not only uses split-reads and discordant read pairs for SV prediction, but also integrates discordant and concordant read pairs (fragments) to genotype SVs under a statistical framework. Specifically, BioSV also provides a multiple-sample-based SV caller for family or population based WGS studies. Moreover, BioSV exhibits high performance on both simulated and real WGS data in SV calling and genotyping.

### <a name="req"></a>Requirements and installation
- Requirements
  
  * *awk*
  
  * *R programming*
  
  * python 2.7 or higher version
  
  * [pysam](https://github.com/pysam-developers/pysam)
  
  * [data.table](https://cran.r-project.org/web/packages/data.table/)
  
  * tools (R package)
  
  * [igraph](https://github.com/igraph/igraph)
  
  * [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

- Download BioSV

`git clone https://github.com/ZhangLiECNU/BioSV.git`

- To install BioSV, users should change directory to "./BioSV", and run 
`sh Install.sh`

Then, all dependencies are prepared and configured, and the executable file 'BioSV' is created.

### <a name="use"></a>Usage 

  Users run the command line `BioSV -h` or `BioSV --help` will show the following messages.
  
  BioSV: Breakpoint-based identification of Structure Variation
            
  Usage: 
`BioSV.main.R  [call/genotype]  [options]`
 
Options:
<pre>
  -b  [file] , --bam= [file] ,  [required]
        
Input bam files (sorted bam files by bwa mem, seperated by comma)
        
  --bedpe [file], [optional]
        
SVs prepared to be genotyped when the mode is set to genotype
        
  -o [directory], --output=[directory]

Directory for output files [BioSV-output]
        
  -t [int], --thread=[int]

Number of threads
        
  -e [flag]

Exclude high coverage regions
        
  --minMRC [int]

Minimal mutant read counts [default = 4]
        
  -h, --help

Print this help message and exit
</pre>

### <a name="exm"></a>Example

- Data preprocessing:
 ```
 bwa mem ref.fa test_1.fq test_2.fq | samtools view -bS - > test.bam
 samtools sort -o test.sort.bam test.bam
 samtools index test.sort.bam
 ```

- Running BioSV for SV calling and genotyping
 ```
 BioSV call -b test.sort.bam -o test.BioSV -t 4 --minMRC 4 
 ```

- Running BioSV only for SV genotyping
 ```
 BioSV genotype -b test.sort.bam --bedpe test.bedpe -o test.BioSV -t 4 --minMRC 4 
 ```

Note: test.bedpe file should be formated as below:

chrom-1 | start-1 | end-1 | chrom-2 | start-2 | end-2 | SV-type | mutant read count (sample-1) | mutant read count (sample-2) | ... | mutant read count (sample-n)
----------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | --------------- | ---------
chr1 |	100	| 101 | chr2 | 202 | 203 | TRA | 20 | 10 |  ...	| 0
chr1 |	200 | 201 | chr1 | 302 | 303 | DEL | 0 | 10 | ... | 20

The order of the samples should be consistent with that of bam files. 

### <a name="tes"></a>Test data

The test data can be simulated by BioSV_simulator. BioSV_simulator is a python script, which can simulate structural variations from diploid genomes, and generate fastq files.
To test the performance of BioSV, users can run test.sh hg19.fasta. The command line will simulate structural variations from chromosomes of 21 and 22, align the simulated reads, and call SV using BioSV. 
If users need to specify any other chromosomes, you can add the names of reference sequences to the file 

`./test_data/chromosomes.txt`

When the script test.sh runs successfully, the BioSV's sensitivity, precision, and genotyping accuracy are recorded in the file of 

`./test_data/performance.txt`



