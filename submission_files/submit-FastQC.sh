# submission script for run-FastQC.sh
# Last updated: 2024-05-23

# NOTES: 
# - argument "-t" needs te be specified as the total number of files to be processed
# - argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-32
#$ -j y
#$ -N QC_before_trim
#$ -m be
#$ -M ummz-arc-records@outlook.com

/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-FastQC.sh /nobackup/ummz/data/transcriptomics/bulk_GCA/run_139_Nov23/renamed /nobackup/ummz/results/transcriptomics/bulk_GCA/all_files_analysis/FastQC/run_139_Nov23 ${SGE_TASK_ID}

#################################################################
# USAGE: 
# (0) /path/to/running/script/run-FastQC.sh 
# (1) /path/to/data/folder 
# (2) /path/to/results/folder 
# (3) ${SGE_TASK_ID}
