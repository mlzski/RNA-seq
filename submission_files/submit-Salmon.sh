# submission script for run-Salmon.sh
# Michal Zulcinski 2022-06-22

# NOTES:
# - argument "-t" needs te be specified as the total number of files to be processed / 2 (both reads are used simultaneously)
# - argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=08:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-80
#$ -j y
#$ -N Salmon_quantif_coh_1_and_2
#$ -m be
#$ -M ummz-arc-records@outlook.com

# define arguments to be passed
script=/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-Salmon.sh
infile=/nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/Salmon_quantification_trimmed/concat_c1_c2/input_files
outfile=//nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/Salmon_quantification_trimmed/concat_c1_c2/
index=/nobackup/ummz/NEW/transcriptomics/results/bulk_GCA/Salmon_quantification_trimmed/concat_c1_c2/salmon_index

# running on transcript level
$script $infile $outfile $index ${SGE_TASK_ID}

###################################################################
# (0) /home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-Salmon.sh
# (1) /path/to/input-data 
# (2) /path/to/output-directory 
# (3) /path/to/index/dir
# (4) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
