# submission script for run-4-picard.sh
# Michal Zulcinski 2020-08-07

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-picard.sh /nobackup/ummz/analyses/run_NEW/4_alignemnt/single-end/bam/outputs /nobackup/ummz/analyses/run_NEW/4_picard/single-end ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/4_picard/single-end/arc_files/output.$JOB_ID.txt
/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-picard.sh /nobackup/ummz/analyses/run_NEW/4_alignemnt/paired-end/bam/outputs /nobackup/ummz/analyses/run_NEW/4_picard/paired-end ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_NEW/4_picard/paired-end/arc_files/output.$JOB_ID.txt

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-picard.sh
# (1) /path/to/data [BAM files]
# (2) /path/to/results 
# (3) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
