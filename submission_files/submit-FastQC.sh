# submission script for run-FastQC.sh
# Last updated: 2023-08-14

# NOTES: 
# - argument "-t" needs te be specified as the total number of files to be processed
# - argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-24
#$ -j y
#$ -N QC_before_trim_concat
#$ -m be
#$ -M ummz-arc-records@outlook.com

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-FastQC.sh /nobackup/ummz/NEW/transcriptomics/data/bulk_GCA/other_data_files/concatenated_samples /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/QC_before_trim/concatenated_samples ${SGE_TASK_ID}

#################################################################
# USAGE: 
# (0) /path/to/running/script/run-FastQC.sh 
# (1) /path/to/data/folder 
# (2) /path/to/results/folder 
# (3) ${SGE_TASK_ID}
