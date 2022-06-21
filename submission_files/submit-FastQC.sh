# submission script for run-FastQC.sh
# Michal Zulcinski 2022-05-25

# NOTES: 
# - argument "-t" needs te be specified as the total number of files to be processed
# - argument "-N" needs to be specified as the name for running job 

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=1G
#$ -t 1-182
#$ -j y
#$ -N QC
#$ -m be
#$ -M ummz-arc-records@outlook.com


# USAGE:
#/path/to/running/script/run-1-qc.sh /path/to/data /path/to/results/ ${SGE_TASK_ID} >> /path/to/arc_files/output.$JOB_ID.txt

/nobackup/ummz/analysis-May-22/1-qc/run-1-qc.sh /nobackup/ummz/analysis-May-22/data /nobackup/ummz/analysis-May-22/1-qc ${SGE_TASK_ID} >> /nobackup/ummz/analysis-May-22/1-qc/arc_files/output.$JOB_ID.txt
