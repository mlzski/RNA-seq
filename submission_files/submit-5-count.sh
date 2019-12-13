# submission script for run-5-count.sh
# Michal Zulcinski 2019-12-12

#$ -cwd -V
#$ -l h_rt=02:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analysis_nov19/RNA-seq/scripts/run-5-count.sh /nobackup/ummz/analysis_nov19/results/4_alignment/bam /nobackup/ummz/analysis_nov19/results/5_counting /nobackup/ummz/reference/annotation/Homo_sapiens.GRCh38.98.gtf ${SGE_TASK_ID} >> /nobackup/ummz/analysis_nov19/results/5_counting/arc_files/output.$JOB_ID.txt

