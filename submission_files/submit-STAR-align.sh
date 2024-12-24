# submission script for run-STAR-align.sh
# Last updated 2023-08-18

# argument "-t" needs te be specified as the total number of samples to be processed [no. of files / 2]
# argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-12
#$ -j y
#$ -N alignment_PE_concat
#$ -m be
#$ -M ummz-arc-records@outlook.com

# for single-end
#/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-STAR-align.sh 'SE' /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/trimming/single_end/concatenated_samples/v1 /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/STAR_alignment_trimmed/alignment/single_end/concatenated_samples/v1 /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/STAR_alignment_trimmed/indexing/index ${SGE_TASK_ID}
 
# for paired-end
/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-STAR-align.sh 'PE' /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/trimming/paired_end/concatenated_samples/v1 /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/STAR_alignment_trimmed/alignment/paired_end/concatenated_samples/v1 /nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/STAR_alignment_trimmed/indexing/index ${SGE_TASK_ID}

#################################################################
# USAGE:
# (0) /home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-STAR-align.sh
# (1) [either SE or PE]
# (2) /path/to/data 
# (3) /path/to/results 
# (4) /path/to/index
# (5) ${SGE_TASK_ID} 
