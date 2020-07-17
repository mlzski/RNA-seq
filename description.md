This file contains a detailed desctiption of the pipeline.

## Project structure:

All permanent files (i.e scripts, tools, readmes, etc.) should be stored in my home directory (/home/home02/ummz).  
All temporary files (i.e. data, reference geome and results) should be placed on /nobackup/ummz, but cannot be stored for a long time. They should be created either before or during the analysis and then removed right after, once they are not needed anymore.  

```
# directory with software
/home/home02/ummz/tools

# directory with reference genome (TO BE CREATED FOR EACH ANALYSIS)
/nobackup/ummz/reference_genome
/home/home02/ummz/reference/README.txt      => this file contains all info about reference genome

# directory with data (TO BE CREATED FOR EACH ANALYSIS; example below)
/nobackup/ummz/analysis_nov19/data

# directory for results (TO BE CREATED FOR EACH ANALYSIS; example below)
/nobackup/ummz/analysis_nov19/results

/home/home02/ummz/RNA-seq  
    description.md  
    LICENSE
    README.md 
    structure.txt
    
    /postprocessing
        postpro-1-qc.sh         => probably obsolete (ckeck it!)
        postpro_qc.py           => probably obsolete (it should be done using ngsReports package in R)
       
    /scripts
        run-1-qc.sh
        run-2-trim.sh
        run-3-qc-trim.sh
        run-4-align.sh
        run-4-index.sh
        run-5-count.sh
        /run_locally             => probably not needed; TO BE REMOVED
    
    /submission_files
        submit-1-qc.sh
        submit-2-trim.sh
        submit-3-qc-trim.sh
        submit-4-align.sh
        submit-4-index.sh
        submit-5-count.sh
```
    
## Notice about reference genome

All reference genome files sholud be downloaded and placed on /nobackup every time before the analysis and then removed right after. It is because we do not want to store them permanently.

## Analysis steps - parameters used in each software (examples):

1) Quality control on raw fastq files (using FastQC)
```
fastqc 
-o /path/to/report/
--threads 4 
--dir /path/to/temp/ 
/path/to/data/[file_name.fastq]
```

2) Trimming (using Trimmomatic)
```
java -jar /home/home02/ummz/tools/Trimmomatic-0.39/trimmomatic-0.39.jar 
PE 
-threads 4 
-phred33 
$read1 
$read2 
/path/to/results/processed_fastq/${coreFile}R1_paired.fq 
/path/to/results/processed_fastq/${coreFile}R1_unpaired.fq 
/path/to/results/processed_fastq/${coreFile}R2_paired.fq 
/path/to/results/processed_fastq/${coreFile}R2_unpaired.fq 
ILLUMINACLIP:/home/home02/ummz/tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 
LEADING:3 
TRAILING:3 
SLIDINGWINDOW:4:15 
MINLEN:36
```

3) Quality control on trimmed fastq files (using FastQC)
```
fastqc 
-o /path/to/report/
--threads 4 
--dir /path/to/temp/ 
/path/to/data/[file_name.fastq]
```

4) Read alignment (using STAR)

NOTICE: it requires 2 new directories to be created and genome (.gtf) and annotation (.fa) files to be downloaded, plus 3 other directories for read mapping: index, sam, bam

```
STAR 
--runThreadN 8 
--runMode genomeGenerate 
--genomeDir ./index 
--genomeFastaFiles ./genome/human_chr22.fa 
--sjdbGTFfile ./annotation/human_chr22.gtf 
--sjdbOverhang 150
```
**NOTICE:** need to unzip fq.gz beforehand

```
STAR 
--runThreadN 8 
--runMode alignReads
--genomeDir [/working/directory]/index 
--readFilesIn [/data/after/trimming]/[sample]_R1_001_val_1.fq.gz [/data/after/trimming]/[sample]_R2_001_val_2.fq.gz 
--outFileNamePrefix [/working/directory]/bam/[sample] 
--outSAMtype BAM SortedByCoordinate 
--outSAMattributes All 
--outSAMstrandField intronMotif
```

# run again to produce the SAM format output
```
STAR 
--runThreadN 8 
--runMode alignReads 
--genomeDir [/working/directory]/index 
--readDilesIn [/data/after/trimming]/[sample]_R1_001_val_1.fq.gz [/data/after/trimming]/[sample]_R2_001_val_2.fq.gz 
--outFileNamePrefix [/working/directory]/sam/[sample] 
--outSAMtype SAM 
--outSAMattributes All 
--outSAMstrandField intronMotif 
```

5) BAM manipulation (using samtools)

indexing using **run-4-samtools.sh**

```
bash run-4-samtools.sh /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_SE/bam /nobackup/ummz/analyses/rerun_FINAL/run_1/samtools_SE.log
bash run-4-samtools.sh /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_PE/bam /nobackup/ummz/analyses/rerun_FINAL/run_1/samtools_PE.log

bash run-4-samtools.sh /nobackup/ummz/analyses/rerun_FINAL/run_2/alignment_SE/bam /nobackup/ummz/analyses/rerun_FINAL/run_2/samtools_SE.log
bash run-4-samtools.sh /nobackup/ummz/analyses/rerun_FINAL/run_2/alignment_PE/bam /nobackup/ummz/analyses/rerun_FINAL/run_2/samtools_PE.log

```

6) featureCounts() 

https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts
