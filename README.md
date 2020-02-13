# RNA-seq

## Analysis steps:

1) Quality control on raw fastq files (using FastQC)
2) Trimming (using Trimmomatic)
3) Quality control on trimmed fastq files (using FastQC)
4) Read alignment (using STAR)
5) BAM manipulation (using samtools) [TO BE IMPLEMENTED]
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
chmod 755 fastqc      # execution rights needs to be given
```

- Trimmomatic (version 0.39) [source: http://www.usadellab.org/cms/?page=trimmomatic]

```
wget "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip"
unzip Trimmomatic-0.39.zip
```

- STAR (version 2.7.3a) [source: https://github.com/alexdobin/STAR/]

```
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
cd STAR-2.7.3a
cd STAR/source        # needs to be compiled
make STAR
```

- samtools (version ) [source:]

```
TBC
```

- Cufflinks (version 2.2.1) [source: http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.OSX_x86_64.tar.gz]

```
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar zxvf cufflinks-2.2.1.tar.gz
```

## Getting started

The pipeline is not automated yet. Each step needs to be launched manually. 

0. Setup: create 5 new directories (one for each step) and copy in running and submission scripts from /home/home02/ummz/RNA-seq/ (e.g run-1-qc.sh & submit-1-qc.sh for the first step)

1. Quality control (FastQC): 
  * create 2 new directories: *report* (for results) and *temp* (for intermediate files) 
  * modify the last line in the submission file: 
```
/path/to/running/script/run-1-qc.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
  * launch the submission file: `qsub submit-[file_name].sh`
  
2. Trimming (Trimmomatic):
 * create a new directory for results *processed_fastq* and 2 sub-directories within it: *paired* and *unpaired*
 * modify the last line in the submission file:
```
/path/to/running/script/run-2-trim.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
  * launch the submission file: `qsub submit-[file_name].sh`
  * **NOTICE:** output files need to be moved to *paired* and *unpaired* folders manually once the running is finished
  
3. Quality control after trimming (FastQC):



4. Alignment (STAR):




5. Reads quantification (Cufflinks):




## Running example 
All the steps were ran on a set of 41 samples.

```
/1_quality_control          => 145 M
/2_trimming                 => 381 G
/3_quality_control_trimmed  => 422 M
/4_alignment                => 111 G
/5_counting                 =>  15 G
```

Files to be backuped:

```
=> all arc_files (from each directory)
=> /1_quality_control/report
=> /3_quality_control_trimmed/report
=> /5_counting
```
