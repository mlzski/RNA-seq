LAST MODIFIED: 13/02/2020

There are 5 steps of the analysis: quality control, trimming, quality control after trimming, alignment and reads quantification. Each step is enclosed in one folder which includes all output files, arc records and 2 shell scripts ("run" and "submit"). The shell scripts are copied from the main directory: /home/home02/ummz/RNA-seq.

**NOTICE 1:** each directory below includes a folder named *arc_files* where all arc output files are stored. It looks similar for each step:
```
  /arc_files
    output.JOB_ID.txt
    submit-[FILE_NAME].sh.eJOB_ID.N (from 1 to N; N => number of files) 
    submit-[FILE_NAME].sh.oJOB_ID.N (from 1 to N; N => number of files)
 ``` 

**NOTICE 2:** there is an estimated size of each folder in braces (e.g. {145 M}) 

**STEP 1: Quality control using FastQC**
```
/1_quality_control {145 M} 

  /arc_files

  /report [1 .html and 1 .zip for each processed file]
    11026_S12_L005_R1_001_fastqc.html (EXAMPLE)
    11026_S12_L005_R1_001_fastqc.zip  (EXAMPLE)

  run-1-qc.sh
  submit-1-qc.sh

  /temp [empty directory]
```

**STEP 2: Trimming using Trimmomatic**

**NOTICE:** */paired* and */unpaired* directories need to be created manualy after trimming and then the files can be moved respectively
```
/2_trimming {381 G}

  /arc_files 
    
  /processed_fastq
    /paired
      11026_S12_L005_R1_paired.fq (EXAMPLE)
      11026_S12_L005_R2_paired.fq (EXAMPLE)
      
    /unpaired
      11026_S12_L005_R1_unpaired.fq (EXAMPLE)
      11026_S12_L005_R2_unpaired.fq (EXAMPLE)
   
  run-2-trim.sh
  submit-2-trim.sh
```

**NOTICE:** there are 41 output files when running Trimmomatics in **SE** mode and 146 when running as **PE**


**STEP 3: Quality control after trimming using FastQC**
```
/3_quality_control_trimmed {422 M}

  /arc_files
  
  /report
    /paired [1 .html and 1 .zip for each processed file]
    /unpaired [1 .html and 1 .zip for each processed file]
  
  run-3-qc-trim.sh
  submit-3-qc-trim.sh

  /temp
```

**STEP 4: Reads alignment using STAR**
```
/4_alignment {111 G}

  /arc_files
    
  /bam [5 output files for each processed file]
    11026_S12_L005_Aligned.sortedByCoord.out.bam  (EXAMPLE) binary file, aligment output sorted by coordinate
    11026_S12_L005_Log.final.out                  (EXAMPLE) summary mapping statistics after mapping job is complete
    11026_S12_L005_Log.out                        (EXAMPLE) log file with information about the run of alignment
    11026_S12_L005_Log.progress.out               (EXAMPLE) contains job progress statistics 
    11026_S12_L005_SJ.out.tab                     (EXAMPLE) contains all counts for each chromosome
  
  Log.out => log file with information about the run of index generation
  run-4-align.sh
  run-4-index.sh      (NOTICE: to be launched once only for a new reference genome)
  submit-4-align.sh
  submit-4-index.sh   (NOTICE: to be launched once only for a new reference genome)
```

**STEP 5: Reads quantification using Cufflinks**

**NOTICE:** the cuffmerge command was used to merge together several Cufflinks assemblies

```
/5_counting {15 G}

  /arc_files
    
  /out_11026_S12_L005 (EXAMPLE) [1 directoriy for each sample (a pair of 2 files)]
    genes.fpkm_tracking     (EXAMPLE) generic FPKM Tracking Format
    isoforms.fpkm_tracking  (EXAMPLE) contains the estimated isoform-level expression values in the generic FPKM Tracking Format
    skipped.gtf             (EXAMPLE) ???
    transcripts.gtf         (EXAMPLE) contains Cufflinks’ assembled isoforms
    
  /results_all_in_one (merged using cuffmerge)
    genes.fpkm_tracking     (EXAMPLE)
    isoforms.fpkm_tracking  (EXAMPLE)
    skipped.gtf             (EXAMPLE)
    transcripts.gtf         (EXAMPLE)
   
  run-5-count.sh
  submit-5-count.sh
```
  

