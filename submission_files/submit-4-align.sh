# submission script for run-4-align.sh
# Michal Zulcinski 2019-12-12

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/path/to/running/script/run-4-align.sh /path/to/input/files /path/to/results /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt

# for single-end
# /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/single-end/run-4-align_SE.sh /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/single-end/processed_fastq /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/single-end /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/single-end/arc_files/output.$JOB_ID.txt

# for paired-end | paired
# /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/paired-end/paired/run-4-align_PE_paired.sh /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/paired-end/processed_fastq/paired /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/paired-end/paired /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/paired-end/paired/arc_files/output.$JOB_ID.txt

# for paired-end | unpaired
# /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/paired-end/unpaired/run-4-align_PE_unpaired.sh /nobackup/ummz/analyses/run_IV_Feb20/2_trimming/paired-end/processed_fastq/unpaired /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/paired-end/unpaired /nobackup/ummz/reference/index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_IV_Feb20/4_alignment/paired-end/unpaired/arc_files/output.$JOB_ID.txt
