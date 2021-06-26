# submission script for run-6-Salmon.sh
# Michal Zulcinski 2021-03-01

#$ -cwd -V
#$ -l h_rt=08:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -t 1-41
#$ -m be
#$ -M ummz-arc-records@outlook.com

#/home/home02/ummz/github_dirs/RNA-seq/scripts/run-6-Salmon_Jun21.sh 'transcript-level' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_PE/processed_fastq/paired/ /nobackup/ummz/analyses/run_17_Jun21/quants_all/transcript-level /nobackup/ummz/analyses/run_17_Jun21/index_all/Homo_sapiens.GRCh37.75_quasi_index /nobackup/ummz/analyses/run_17_Jun21/index_all/Homo_sapiens.GRCh37.75.gtf ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_17_Jun21/arc_files/transcript-level/output.$JOB_ID.txt

/home/home02/ummz/github_dirs/RNA-seq/scripts/run-6-Salmon_Jun21.sh 'gene-level' /nobackup/ummz/analyses/run_12_Aug20/2_trimming_PE/processed_fastq/paired/ /nobackup/ummz/analyses/run_17_Jun21/quants_all/gene-level /nobackup/ummz/analyses/run_17_Jun21/index_all/Homo_sapiens.GRCh37.75_quasi_index /nobackup/ummz/analyses/run_17_Jun21/index_all/Homo_sapiens.GRCh37.75.gtf ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_17_Jun21/arc_files/gene-level/output.$JOB_ID.txt

 
###################################################################
# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-6-Salmon_Jun21.sh
# (1) [transcript-level or gene-level]
# (2) /path/to/input-data 
# (3) /path/to/output-directory 
# (4) /path/to/index/dir
# (5) /path/to/gtf_file (only used for 'gene-level')
# (6) ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
