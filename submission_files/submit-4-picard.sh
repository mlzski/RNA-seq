# submission script for run-4-picard.sh
# Michal Zulcinski 2020-08-07

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-4-picard.sh /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_all/bam /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_all ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_all/arc_files/output.$JOB_ID.txt

#/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-4-picard.sh /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_noXY/bam /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_noXY ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_noXY/arc_files/output.$JOB_ID.txt

#/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-4-picard.sh /nobackup/ummz/analyses/run_12_Aug20/4_alignment_PE_all/bam /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_all ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_all/arc_files/output.$JOB_ID.txt

#/home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-4-picard.sh /nobackup/ummz/analyses/run_12_Aug20/4_alignment_PE_noXY/bam /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_noXY ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_noXY/arc_files/output.$JOB_ID.txt


# (0) /home/home02/ummz/GitHubRepos/RNA-seq/scripts/run-4-picard.sh
# (1) /path/to/data [BAM files]
# (2) /path/to/results 
# (3) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
