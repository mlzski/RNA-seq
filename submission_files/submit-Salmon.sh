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

# define arguments to be passed
script=/nobackup/ummz/analysis-May-22/read-quantification/run-Salmon.sh
infile=/nobackup/ummz/analysis-May-22/2-trim-compressed/trimmed/paired-end/paired
outfile=/nobackup/ummz/analysis-May-22/read-quantification/output-quantif
index=/nobackup/ummz/analysis-May-22/reference-genome/salmon_index
gtf=/nobackup/ummz/analysis-May-22/reference-genome/gencode.v40.annotation.gtf.gz

# running on transcript level (gtf file is not used)
$script 'transcript-level' $infile $outfile $index $gtf ${SGE_TASK_ID} >> /nobackup/ummz/analysis-May-22/read-quantification/arc_files/output-transcript-level.$JOB_ID.txt

# running on gene level
### $script 'gene-level' $infile $outfile $index $gtf ${SGE_TASK_ID} >> /nobackup/ummz/analysis-May-22/read-quantification/arc_files/output-gene-level.$JOB_ID.txt

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
