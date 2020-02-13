# submission script for run-5-count.sh
# Michal Zulcinski 2019-12-12

#$ -cwd -V
#$ -l h_rt=02:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-5-count.sh /path/to/files/bam /path/to/results /nobackup/ummz/reference/annotation/Homo_sapiens.GRCh38.98.gtf ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt
