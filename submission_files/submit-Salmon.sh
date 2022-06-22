# submission script for run-Salmon.sh
# Michal Zulcinski 2022-06-22

# NOTES:
# - argument "-t" needs te be specified as the total number of files to be processed / 2 (both reads are used simultaneously)
# - argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=08:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-91
#$ -j y
#$ -N read-quantif
#$ -m be
#$ -M ummz-arc-records@outlook.com

 
### running on transcript level (gtf file is not used)

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-Salmon.sh \
'transcript-level' \
/nobackup/ummz/analysis-May-22/2-trim-compressed/test \
/nobackup/ummz/analysis-May-22/read-quantification/output \
/nobackup/ummz/analysis-May-22/reference-genome/salmon_index \
/nobackup/ummz/analysis-May-22/reference-genome/gencode.v40.annotation.gtf.gz \
${SGE_TASK_ID} >> /nobackup/ummz/analysis-May-22/read-quantification/arc_files/output-transcript-level.$JOB_ID.txt


### running on gene level

#/home/home02/ummz/github_dirs/RNA-seq/scripts/run-Salmon.sh \
#'gene-level' \
#/nobackup/ummz/analysis-May-22/2-trim-compressed/test \
#/nobackup/ummz/analysis-May-22/read-quantification/output \
#/nobackup/ummz/analysis-May-22/reference-genome/salmon_index \
#/nobackup/ummz/analysis-May-22/reference-genome/gencode.v40.annotation.gtf.gz \
#${SGE_TASK_ID} >> /nobackup/ummz/analysis-May-22/read-quantification/arc_files/output-gene-level.$JOB_ID.txt

###################################################################
# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-Salmon.sh
# (1) [transcript-level or gene-level]
# (2) /path/to/input-data 
# (3) /path/to/output-directory 
# (4) /path/to/index/dir
# (5) /path/to/gtf_file (only used for 'gene-level')
# (6) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
