# submission script for run-4-picard.sh
# Michal Zulcinski 2020-05-17

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analyses/rerun_Ian/rerun_5/scripts/run-4-picard.sh /nobackup/ummz/analyses/rerun_Ian/rerun_5/bam/outputs /nobackup/ummz/analyses/rerun_Ian/rerun_5/pic ${SGE_TASK_ID} >> /nobackup/ummz/analyses/rerun_Ian/rerun_5/output.$JOB_ID.txt

