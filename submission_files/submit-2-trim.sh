# submission script for run-2-trim.sh
# Michal Zulcinski 2019-12-11

#$ -cwd -V
#$ -l h_rt=1:00:00
#$ -l h_vmem=1G
#$ -l np=4
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-2-trim.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt

# for single-end
# /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/single-end/run-2-trim.sh /nobackup/ummz/analyses/data /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/single-end/ ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/single-end/arc_files/output.$JOB_ID.txt

# for paired-end
# /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/paired-end/run-2-trim.sh /nobackup/ummz/analyses/data /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/paired-end/ ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/paired-end/arc_files/output.$JOB_ID.txt
