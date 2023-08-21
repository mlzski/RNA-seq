# submission script for run-STAR-index.sh
# Last updated 2023-08-18

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -j y
#$ -N indexing_STAR
#$ -m be
#$ -M ummz-arc-records@outlook.com

bash run-STAR-index.sh index reference_genome/GRCh38.primary_assembly.genome.fa reference_genome/gencode.v44.primary_assembly.annotation.gtf >> STAR_index_${JOB_ID}.log

#################################################################
# USAGE:
# (0) /path/to/running/script/run-STAR-index.sh 
# (1) /path/to/folder/for/index 
# (2) /path/to/genome.fa 
# (3) /path/to/annotation.gtf 
# (4) >> /path/to/arc_files/output.$JOB_ID.txt
