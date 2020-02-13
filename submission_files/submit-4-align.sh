# submission script for run-4-align.sh
# Michal Zulcinski 2019-12-12

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-4-align.sh /path/to/input/files /path/to/results /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
