# RNA-seq

## Analysis steps:

1) Quality control on raw fastq files (using FastQC)
2) Trimming (using Trimmomatic)
3) Quality control on trimmed fastq files (using FastQC)
4) Read alignment (using STAR)
5) BAM manipulation (using samtools)
6) Obtain read counts (using Cufflinks)

## Download

Clone this repository to your own directory

```
git clone https://github.com/mihaux/RNA-seq.git
```

## Requirements

1) Access to a Linux-based HPC service. (Tested on ARC4 cluster, based on the CentOS 7 distribution)

2) Reference genome. It can be downloaded from https://genome.ucsc.edu/index.html

3) Following tools need to be downloaded and installed prior to running the pipeline:

- FastQC (version 0.11.8) [source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/]

```
wget "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip"

unzip fastqc_v0.11.8.zip

cd FastQC

chmod 755 fastqc
```
- Trimmomatic (version 0.39) [source: http://www.usadellab.org/cms/?page=trimmomatic]

```
wget "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"

unzip Trimmomatic-0.39.zip
```
- STAR (version ) [source:]

- samtools (version ) [source:]

- Cufflinks () [source:]


## Getting started
