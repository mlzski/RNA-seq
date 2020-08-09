# submission script for run-5-counting.R
# Michal Zulcinski 2020-03-06

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -t 1-2
#$ -m be
#$ -M ummz-arc-records@outlook.com

# Rscript /home/home02/ummz/RNA-seq/scripts/run-5-featureCounts_I.R /nobackup/ummz/analyses/run_I_Nov19/4_aligment/bam /nobackup/ummz/analyses/run_I_Nov19/5_featureCounts/postprocessed ${SGE_TASK_ID} >> /nobackup/ummz/analyses/run_I_Nov19/5_featureCounts/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/RNA-seq/scripts/run-5-featureCounts_Ian.R dups SE /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_1/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_1/arc_files/SE/output.$JOB_ID.txt
#Rscript /home/home02/ummz/RNA-seq/scripts/run-5-featureCounts_Ian.R nodups SE /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_1/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_1/arc_files/SE/output.$JOB_ID.txt

#Rscript /home/home02/ummz/RNA-seq/scripts/run-5-featureCounts_Ian.R dups SE /nobackup/ummz/analyses/rerun_FINAL/run_2/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_2/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_2/arc_files/SE/output.$JOB_ID.txt
#Rscript /home/home02/ummz/RNA-seq/scripts/run-5-featureCounts_Ian.R nodups SE /nobackup/ummz/analyses/rerun_FINAL/run_2/alignment_SE/bam/outputs /nobackup/ummz/reference_genome/NEW/annotation/hg38.refGene.gtf /nobackup/ummz/analyses/rerun_FINAL/run_2/featCounts_SE >> /nobackup/ummz/analyses/rerun_FINAL/run_2/arc_files/SE/output.$JOB_ID.txt

# this line needs to be modified every time:  #$ -t 1-2
# this line was different in Ian's script:    #$ -l h_rt=04:00:00

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.sh
# (1) ??? /path/to/data 
# (2) ??? /path/to/results 
# (3) ??? ${SGE_TASK_ID} 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt
