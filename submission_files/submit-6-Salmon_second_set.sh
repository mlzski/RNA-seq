# submission script for run-6-Salmon.sh
# Michal Zulcinski 2021-03-01

#$ -cwd -V
#$ -l h_rt=08:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-12
#$ -m be
#$ -M ummz-arc-records@outlook.com

# single-end [TO BE MODIFIED]
#/home/home02/ummz/github_dirs/RNA-seq/scripts/run-4-align.sh 'SE' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_SE/processed_fastq /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_all/ /nobackup/ummz/reference_genome/generatedBySTAR_July20/index_allChr ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_12_Aug20/4_alignment_SE_all/arc_files/output.$JOB_ID.txt


# paired-end
/home/home02/ummz/github_dirs/RNA-seq/scripts/run-6-Salmon.sh 'PE' /nobackup/ummz/analyses/run_15_Mar21/second_set_INPUT/ /nobackup/ummz/analyses/run_15_Mar21/quants_second_set /nobackup/ummz/analyses/run_13_Jan21/index/Homo_sapiens.GRCh38_index ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_15_Mar21/arc_files/output.$JOB_ID.txt

 
###################################################################
# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-6-Salmon.sh
# (1) [SE or PE]
# (2) /path/to/input-data 
# (3) /path/to/output-directory 
# (4) /path/to/index
# (5) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
