# submission script for run-3-qc-trim.sh
# Michal Zulcinski 2019-12-11

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-164
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-3-qc-trim.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt

