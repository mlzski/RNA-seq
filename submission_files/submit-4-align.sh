# submission script for run-4-align.sh
# Michal Zulcinski 2020-08-07

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh 'SE' /nobackup/ummz/analyses/run_NEW/2_trimming/single-end/processed_fastq /nobackup/ummz/analyses/run_NEW/4_alignment/single-end /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/4_alignment/single-end/arc_files/output.$JOB_ID.txt
/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh 'PE' /nobackup/ummz/analyses/run_NEW/2_trimming/paired-end/processed_fastq /nobackup/ummz/analyses/run_NEW/4_alignment/paired-end /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/4_alignment/paired-end/arc_files/output.$JOB_ID.txt

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh
# (1) [either SE or PE]
# (2) /path/to/data 
# (3) /path/to/results 
# (4) /path/to/index
# (5) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
