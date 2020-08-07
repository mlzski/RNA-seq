# submission script for run-2-trim.sh
# Michal Zulcinski 2020-08-07

#$ -cwd -V
#$ -l h_rt=1:00:00
#$ -l h_vmem=1G
#$ -l np=4
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-2-trim.sh PE /nobackup/ummz/analyses/data /nobackup/ummz/analyses/run_NEW/2_trimming/paired-end/ ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/2_trimming/paired-end/arc_files/output.$JOB_ID.txt

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-2-trim.sh
# (1) [either SE or PE]
# (2) /path/to/data 
# (3) /path/to/results 
# (4) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
