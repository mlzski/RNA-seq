# submission script for run-3-qc-trim.sh
# Michal Zulcinski 2020-08-07

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-84.  # TO BE CHECKED (EITER 82 OR 164)
#$ -m be
#$ -M ummz-arc-records@outlook.com

# for single-end
# /home/home02/ummz/github_dirs/RNA-seq/scripts/run-3-qc-trim.sh /nobackup/ummz/analyses/run_NEW/2_trimming/single-end/processed_fastq/ /nobackup/ummz/analyses/run_NEW/3_quality_control_trimmed/single-end ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/3_quality_control_trimmed/single-end/arc_files/output.$JOB_ID.txt

# for paired-end (no need to run for unpaired)
# /home/home02/ummz/github_dirs/RNA-seq/scripts/run-3-qc-trim.sh /nobackup/ummz/analyses/run_NEW/2_trimming/paired-end/processed_fastq/ /nobackup/ummz/analyses/run_NEW/3_quality_control_trimmed/paired-end ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/3_quality_control_trimmed/paired-end/arc_files/output.$JOB_ID.txt

# NOTICE: no need to launch for 'unpaired' files

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-3-qc-trim.sh
# (1) /path/to/data 
# (2) /path/to/results 
# (3) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
