# submission script for run-5-featureCounts_Ian.R
# Michal Zulcinski 2020-05-20

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -m be
#$ -M ummz-arc-records@outlook.com

#Rscript /nobackup/ummz/analyses/rerun_FINAL/run_1/run-5-featureCounts_Ian.R dups SE /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_1/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_1/arc_files/SE/output.$JOB_ID.txt
#Rscript /nobackup/ummz/analyses/rerun_FINAL/run_1/run-5-featureCounts_Ian.R nodups SE /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_1/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_1/arc_files/SE/output.$JOB_ID.txt

#Rscript /nobackup/ummz/analyses/rerun_FINAL/run_1/run-5-featureCounts_Ian.R dups SE /nobackup/ummz/analyses/rerun_FINAL/run_2/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_2/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_2/arc_files/SE/output.$JOB_ID.txt
#Rscript /nobackup/ummz/analyses/rerun_FINAL/run_1/run-5-featureCounts_Ian.R nodups SE /nobackup/ummz/analyses/rerun_FINAL/run_2/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_2/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_2/arc_files/SE/output.$JOB_ID.txt
