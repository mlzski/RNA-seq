# submission script for run-4-align.sh
# Michal Zulcinski 2019-12-11

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-164
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analysis_nov19/RNA-seq/scripts/run-3-qc-trim.sh /nobackup/ummz/analysis_nov19/results/2_trimming/processed_fastq /nobackup/ummz/analysis_nov19/results/3_quality_control_trimmed ${SGE_TASK_ID} >> /nobackup/ummz/analysis_nov19/results/3_quality_control_trimmed/arc_files/output.$JOB_ID.txt


