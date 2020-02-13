# submission script for run-1-qc.sh
# Michal Zulcinski 2019-12-05

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-82
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-1-qc.sh /path/to/data /path/to/results/ ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
