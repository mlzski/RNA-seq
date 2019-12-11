# submission script for run-1-qc.sh
# Michal Zulcinski 2019-12-05

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-82
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analysis_nov19/RNA-seq/scripts/run-1-qc.sh /nobackup/ummz/analysis_nov19/data /nobackup/ummz/analysis_nov19/results/1_quality_control ${SGE_TASK_ID} >> /nobackup/ummz/analysis_nov19/results/1_quality_control/arc_files/output.$JOB_ID.txt

