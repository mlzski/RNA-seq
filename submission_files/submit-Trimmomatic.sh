# submission script for run-Trimmomatic.sh
# Michal Zulcinski 2020-06-07

# NOTES:
# - argument "-t" needs te be specified as the total number of files to be processed / 2 (both reads are used simultaneously)
# - argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=1:00:00
#$ -l h_vmem=1G
#$ -l np=4
#$ -t 1-91
#$ -j y
#$ -N trimming
#$ -m be
#$ -M ummz-arc-records@outlook.com


# USAGE:
# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-2-trim.sh
# (1) [either 'SE' or 'PE']
# (2) /path/to/data 
# (3) /path/to/results 
# (4) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt

/nobackup/ummz/analysis-May-22/2-trim/run-2-trim.sh 'PE' /nobackup/ummz/analysis-May-22/data /nobackup/ummz/analysis-May-22/2-trim ${SGE_TASK_ID} >> /nobackup/ummz/analysis-May-22/2-trim/arc_files/output.$JOB_ID.txt

