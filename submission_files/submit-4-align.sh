# submission script for run-4-align.sh
# Michal Zulcinski 2019-12-12

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analysis_nov19/RNA-seq/scripts/run-4-align.sh /nobackup/ummz/analysis_nov19/results/2_trimming/processed_fastq/paired /nobackup/ummz/analysis_nov19/results/4_alignment /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /nobackup/ummz/analysis_nov19/results/4_alignment/arc_files/output.$JOB_ID.txt

