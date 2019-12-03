This file contains a detailed desctiption of the pipeline.


## Project structure:

```
# directory with software
/nobackup/ummz/tools

# directory with reference genome
/nobackup/ummz/reference_genome

# directory with data
/nobackup/ummz/analysis_nov19/data

~/RNA-seq  
    description.md  
    LICENSE  
    README.md  
    /results  
    /scripts  
    /submission_files  
```
    
    
    
## Analysis steps:


1) Quality control on raw fastq files (using FastQC) => script: run_qc.sh

```
export PATH=/nobackup/leedsomics_workshop/tools/FastQC:$PATH

fastqc -o /nobackup/ummz/analysis_nov19/1_quality_control/report/ --threads 4 --dir /nobackup/ummz/analysis_nov19/1_quality_control/temp/ /nobackup/ummz/analysis_nov19/data/data/[file_name.fastq]
```

NOTICE: it requires two folders: /report for .html report files and /temp for temporary files which are removed before the process ends (so no output in this folder)

2) Trimming (using Trimmomatic)




3) Quality control on trimmed fastq files (using FastQC)

4) Read alignment (using STAR)

```
export PATH=/nobackup/leedsomics_workshop/tools/STAR-2.7.0a/bin/Linux_x86_64_static:$PATH

# it requires to create 2 new directorie and download genome (.gtf) and annotation (.fa) files
# plus 3 other directories for read mapping: index, sam, bam

# current directory: /nobackup/ummz/analysis_nov19/4_read_alignment

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ./genome/human_chr22.fa --sjdbGTFfile ./annotation/human_chr22.gtf --sjdbOverhang 150

# need to unzip fq.gz beforehand

STAR --runThreadN 4 --runMode alignReads --genomeDir [/working/directory]/index --readFilesIn [/data/directory/after/trimming]/[sample]_R1_001_val_1.fq.gz [/data/directory/after/trimming]/[sample]_R2_001_val_2.fq.gz --outFileNamePrefix [/working/directory]/bam/[sample] --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif

# run again to produce the SAM format output

STAR --runThreadN 4 --runMode alignReads --genomeDir [/working/directory]/index --readDilesIn [/data/directory/after/trimming]/[sample]_R1_001_val_1.fq.gz [/data/directory/after/trimming]/[sample]_R2_001_val_2.fq.gz --outFileNamePrefix [/working/directory]/sam/[sample] --outSAMtype SAM --outSAMattributes All --outSAMstrandField intronMotif 

```

5) BAM manipulation (using samtools)



6) Obtain read counts (using Cufflinks)



