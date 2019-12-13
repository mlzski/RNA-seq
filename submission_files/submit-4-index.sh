# submission script for run-4-index.sh
# Michal Zulcinski 2019-12-12

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=9G
#$ -pe smp 8
#$ -m be
#$ -M ummz-arc-records@outlook.com

/nobackup/ummz/analysis_nov19/RNA-seq/scripts/run-4-index.sh /nobackup/ummz/reference/index /nobackup/ummz/reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa /nobackup/ummz/reference/annotation/Homo_sapiens.GRCh38.98.gtf >> /nobackup/ummz/analysis_nov19/results/4_alignment/arc_files/output.$JOB_ID.txt


