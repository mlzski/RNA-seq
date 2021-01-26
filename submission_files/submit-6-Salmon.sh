# submission script for run-4-align.sh
# Michal Zulcinski 2020-08-07

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

# [TO BE MODIFIED]
#/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh 'SE' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_SE/processed_fastq /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_all/ /nobackup/ummz/reference_genome/generatedBySTAR_July20/index_allChr ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_all/arc_files/output.$JOB_ID.txt

# [TO BE MODIFIED]
#/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh 'SE' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_SE/processed_fastq /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_noXY/ /nobackup/ummz/reference_genome/generatedBySTAR_July20/index_noChrXandY ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_noXY/arc_files/output.$JOB_ID.txt

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-6-Salmon.R 'PE' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_PE/processed_fastq/paired /nobackup/ummz/analyses/run_13_Jan21/quants /nobackup/ummz/analyses/run_13_Jan21/index/Homo_sapiens.GRCh38_index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_13_Jan21/arc_files/output.$JOB_ID.txt
 
#/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh 'PE' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_PE/processed_fastq/paired /nobackup/ummz/analyses/run_12_Aug20/4_alignment_PE_noXY/ /nobackup/ummz/reference_genome/generatedBySTAR_July20/index_noChrXandY ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_alignment_PE_noXY/arc_files/output.$JOB_ID.txt

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh
# (1) [either SE or PE]
# (2) /path/to/data 
# (3) /path/to/results 
# (4) /path/to/index
# (5) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
