# RNA-seq

## Analysis steps:

1) Quality control on raw fastq files (using FastQC)
2) Trimming (using Trimmomatic)
3) Quality control on trimmed fastq files (using FastQC)
4) Read alignment (using STAR, bowtie2 can be worth trying out)
5) BAM manipulation (using samtools) [TO BE IMPLEMENTED]
6) Obtain read counts (using Cufflinks or featureCounts)

**NOTICE:** Cufflinks turned out to be dedicatedfor transcript discovery, therefore a new software "featureCounts" was proposed

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

**Setup:** create 5 new directories (one for each step) and copy in running and submission scripts from /home/home02/ummz/RNA-seq/ (e.g run-1-qc.sh & submit-1-qc.sh for the first step)

**1. Quality control (FastQC):**  
  * modify the last line in the submission file: 
```
/path/to/running/script/run-1-qc.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
  * launch the submission file: `qsub submit-1-qc.sh`
  
**2. Trimming (Trimmomatic):**
 * create a new directory for results *processed_fastq* and 2 sub-directories within the *report* folder: *paired* and *unpaired*
 * modify the last line in the submission file:
```
/path/to/running/script/run-2-trim.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
  * launch the submission file: `qsub submit-2-trim.sh`
  * **NOTICE:** output files need to be moved to *paired* and *unpaired* folders manually once the running is finished
  
**3. Quality control after trimming (FastQC):**
  * create 2 new directories: *report* (for results) and *temp* (for intermediate files) and 2 sub-directories within it: *paired* and *unpaired*
  * modify the last line in the submission file: 
```
/path/to/running/script/run-3-qc-trim.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
  * launch the submission file: `qsub submit-3-qc-trim.sh`
  * **NOTICE:** output files need to be moved to *paired* and *unpaired* folders manually once the running is finished

**4. Alignment (STAR):**
 * NOTICE: when running for the first time, follow the instruction from /home/home02/ummz/reference/README.txt (i.e. generate index first) 
 * modify the last line in the INDEX submission file: (first time running only) 
```
/path/to/running/script/run-4-index.sh /nobackup/ummz/reference/index /nobackup/ummz/reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa /nobackup/ummz/reference/annotation/Homo_sapiens.GRCh38.98.gtf >> /path/to/arc_files/output.$JOB_ID.txt
```
 * launch the submission file: `qsub submit-4-index.sh` (first time running only) 
 
 * create a new directory for results *bam*
 * modify the last line in the ALIGN submission file: 
```
/path/to/running/script/run-4-align.sh /path/to/input/files /path/to/results /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
 * launch the submission file: `qsub submit-4-align.sh`

**5. Reads quantification (2 software):**

**5.A. Cufflinks:**
 * modify the last line in the submission file:
 ```
/path/to/running/script/run-5-count.sh /path/to/files/bam /path/to/results /nobackup/ummz/reference/annotation/Homo_sapiens.GRCh38.98.gtf ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
```
 * launch the submission file: `qsub submit-5-count.sh`

**5.B. featureCounts:**
in R

Before launching:
**module load anaconda**
**source activate r36**

**NOTICE:** use 'conda env list' to see all available environments


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
=> all arc_files folders (from each directory)
=> /1_quality_control/report
=> /3_quality_control_trimmed/report
=> /4_alignment/Log.out
=> /4_alignment/[sample_name]_Log.final.out (for all samples)
=> /5_counting (all???)
```
