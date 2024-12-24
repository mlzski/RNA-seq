# submission script for run-Trimmomatic.sh
# Last updated: 2024-05-23

# NOTES:
# - argument "-t" needs te be specified as the total number of samples to be processed [total number of files / 2]
# - argument "-N" needs to be specified as the name for running job

#$ -cwd -V
#$ -l h_rt=2:00:00
#$ -l h_vmem=1G
#$ -t 1-16
#$ -j y
#$ -N trimming_PE_run_139
#$ -m be
#$ -M ummz-arc-records@outlook.com

# for single-end
#/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-Trimmomatic.sh 'SE' /nobackup/ummz/NEW/transcriptomics/data/bulk_GCA/other_data_files/concatenated_samples /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/trimming/single_end/concatenated_samples/v1 ${SGE_TASK_ID}

# for paired-end
/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-Trimmomatic.sh 'PE' /nobackup/ummz/data/transcriptomics/bulk_GCA/run_139_Nov23/renamed/ /nobackup/ummz/results/transcriptomics/bulk_GCA/trimming/paired_end/run_139_Nov23/v1 ${SGE_TASK_ID}

#################################################################
# USAGE:
# (0) /home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-Trimmomatic.sh
# (1) [either 'SE' or 'PE']
# (2) /path/to/input/data 
# (3) /path/to/results 
# (4) ${SGE_TASK_ID} 
