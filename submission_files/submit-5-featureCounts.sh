# submission script for run-5-featureCounts.R
# Michal Zulcinski 2020-03-06

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -t 1-2
#$ -m be
#$ -M ummz-arc-records@outlook.com

Rscript /home/home02/ummz/RNA-seq/scripts/run-5-featureCounts.R /nobackup/ummz/analyses/run_I_Nov19/4_aligment/bam /nobackup/ummz/analyses/run_I_Nov19/5_featureCounts/postprocessed ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_I_Nov19/5_featureCounts/arc_files/output.$JOB_ID.txt

# Rscript /path/to/running/script/run-5-featureCounts.R /path/to/files/bam /path/to/results /nobackup/ummz/reference/annotation/Homo_sapiens.GRCh38.98.gtf ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
