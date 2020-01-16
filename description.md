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

## Analysis steps - parameters used in each software:

1) Quality control on raw fastq files (using FastQC) => scripts: run-1-qc.sh & submit-1-qc.sh

```
fastqc 
-o /nobackup/ummz/analysis_nov19/1_quality_control/report/ 
--threads 4 
--dir /nobackup/ummz/analysis_nov19/1_quality_control/temp/ 
/nobackup/ummz/analysis_nov19/data/data/[file_name.fastq]
```

NOTICE: it requires two folders: /report for .html report files and /temp for temporary files which are removed before the process ends (so no output in this folder)

2) Trimming (using Trimmomatic)

```

```

3) Quality control on trimmed fastq files (using FastQC)

```

```

4) Read alignment (using STAR)

NOTICE: it requires 2 new directories to be created and genome (.gtf) and annotation (.fa) files to be downloaded, plus 3 other directories for read mapping: index, sam, bam

```
STAR 
--runThreadN 4 
--runMode genomeGenerate 
--genomeDir ./index 
--genomeFastaFiles ./genome/human_chr22.fa 
--sjdbGTFfile ./annotation/human_chr22.gtf 
--sjdbOverhang 150

# NOTICE: need to unzip fq.gz beforehand

STAR 
--runThreadN 4 
--runMode alignReads
--genomeDir [/working/directory]/index 
--readFilesIn [/data/after/trimming]/[sample]_R1_001_val_1.fq.gz [/data/after/trimming]/[sample]_R2_001_val_2.fq.gz 
--outFileNamePrefix [/working/directory]/bam/[sample] 
--outSAMtype BAM SortedByCoordinate 
--outSAMattributes All 
--outSAMstrandField intronMotif

# run again to produce the SAM format output

STAR 
--runThreadN 4 
--runMode alignReads 
--genomeDir [/working/directory]/index 
--readDilesIn [/data/after/trimming]/[sample]_R1_001_val_1.fq.gz [/data/after/trimming]/[sample]_R2_001_val_2.fq.gz 
--outFileNamePrefix [/working/directory]/sam/[sample] 
--outSAMtype SAM 
--outSAMattributes All 
--outSAMstrandField intronMotif 
```

5) BAM manipulation (using samtools)



6) Obtain read counts (using Cufflinks)



