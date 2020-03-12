# submission script for run-3-qc-trim.sh
# Michal Zulcinski 2019-12-11

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-164
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-3-qc-trim.sh /path/to/data /path/to/results ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt

# for single-end
#/nobackup/ummz/analyses/run_IV_Feb20/3_quality_control_trimmed/run-3-qc-trim.sh /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/single-end/processed_fastq/ /nobackup/ummz/analyses/run_IV_Feb20/3_quality_control_trimmed/single-end ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/3_quality_control_trimmed/single-end/arc_files/output.$JOB_ID.txt
# for paired-end
#/nobackup/ummz/analyses/run_IV_Feb20/3_quality_control_trimmed/run-3-qc-trim.sh /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/paired-end/processed_fastq/ /nobackup/ummz/analyses/run_IV_Feb20/3_quality_control_trimmed/paired-end ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/3_quality_control_trimmed/paired-end/arc_files/output.$JOB_ID.txt
