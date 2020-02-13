/1_quality_control {145 M}
  /arc_files  [165 files]
    output.24845.txt
    submit-1-qc.sh.e24845.1 (from 1 to 82)
    submit-1-qc.sh.o24845.1 (from 1 to 82)
    
  /postprocessed [to be changed, use ngsReports package in R]
    all_merged.txt          => all outputs files merged together (to be removed/modified)
    qc_metrics_out.txt      => output from postprocessing file (to be removed/modified)
  
    /unzipped [82 files; includes all .zip files from /report]
      11026_S12_L005_R1_001_fastqc 
      11026_S12_L005_R2_001_fastqc

  /report [164 files: 41 .html for each read; 41 .zip for each read]
    11026_S12_L005_R1_001_fastqc.html 
    11026_S12_L005_R1_001_fastqc.zip 
    
  /temp [empty directory]

/2_trimming {381 G}
  /arc_files [83 files]
    output.26787.txt            
    submit-2-trim.sh.e26787.1 (from 1 to 41) 
    submit-2-trim.sh.o26787.1 (from 1 to 41)
    
   /processed_fastq
    /paired [82 files]
      11026_S12_L005_R1_paired.fq
      11026_S12_L005_R2_paired.fq
      
    /unpaired [82 files]
      11026_S12_L005_R1_unpaired.fq
      11026_S12_L005_R2_unpaired.fq

/3_quality_control_trimmed {422 M}
  /arc_files [329 files]
    output.26802.txt
    submit-3-qc-trim.sh.e26802.1 (from 1 to 164)
    submit-3-qc-trim.sh.e26802.1 (from 1 to 164)
    
  /postprocessed
    all_merged_paired.txt
    all_merged_unpaired.txt
    qc_metrics_step3_paired_out.txt
    qc_metrics_step3_unpaired_out.txt
    /unzipped
      /paired [82 directories]
      /unpaired [82 directories]
  
  /report
    /paired [164 files: 41 .html for each read; 41 .zip for each read]
    /unpaired [164 files: 41 .html for each read; 41 .zip for each read]
  
  /temp

/4_alignment {111 G}
  /arc_files
    Log.out
    output.27047.txt
    output.27063.txt
    submit-4-align.sh.e27063.1 (from 1 to 41)
    submit-4-align.sh.o27063.1 (from 1 to 41)
    submit-4-index.sh.e27047
    submit-4-index.sh.o27047
    
  /bam [205 files; 5 output files for each sample]
    11026_S12_L005_Aligned.sortedByCoord.out.bam
    11026_S12_L005_Log.final.out
    11026_S12_L005_Log.out
    11026_S12_L005_Log.progress.out
    11026_S12_L005_SJ.out.tab
    
/5_counting {15 G}
  /arc_files [166 files]
    output.27114.txt
    output.27118.txt
    submit-5-count.sh.e27114.1 (from 1 to 41)
    submit-5-count.sh.o27114.1 (from 1 to 41)
    submit-5-count.sh.e27118.1 (from 1 to 41)
    submit-5-count.sh.e27118.1 (from 1 to 41)
    
  /out_11026_S12_L005 [41 directories like this (for each sample)]
    genes.fpkm_tracking  
    isoforms.fpkm_tracking  
    skipped.gtf  
    transcripts.gtf
    
  /results_all_in_one
    genes.fpkm_tracking
    isoforms.fpkm_tracking
    skipped.gtf
    transcripts.gtf
  
  

