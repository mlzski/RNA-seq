# submission script for run-2-trim.sh
# Michal Zulcinski 2019-12-11

#$ -cwd -V
#$ -l h_rt=1:00:00
#$ -l h_vmem=1G
#$ -l np=4
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analysis_nov19/RNA-seq/scripts/run-2-trim.sh /nobackup/ummz/analysis_nov19/data /nobackup/ummz/analysis_nov19/results/2_trimming ${SGE_TASK_ID} >> /nobackup/ummz/analysis_nov19/results/2_trimming/arc_files/output.$JOB_ID.txt

