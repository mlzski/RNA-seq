# submission script for run-5-counting.R
# Michal Zulcinski 2020-08-10

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -m be
#$ -M ummz-arc-records@outlook.com

Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R dups SE /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_all /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_all >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_all/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R nodups SE /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_all /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_all >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_all/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R dups SE /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_noXY /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_noXY >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_noXY/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R nodups SE /nobackup/ummz/analyses/run_12_Aug20/4_picard_SE_noXY /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_noXY >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_SE_noXY/arc_files/output.$JOB_ID.txt

#---

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R dups PE /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_all /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_all >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_all/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R nodups PE /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_all /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_all >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_all/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R dups PE /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_noXY /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_noXY >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_noXY/arc_files/output.$JOB_ID.txt

#Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.R nodups PE /nobackup/ummz/analyses/run_12_Aug20/4_picard_PE_noXY /nobackup/ummz/reference_genome/generatedBySTAR_July20/annotation/hg38.ncbiRefSeq.gtf /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_noXY >> /nobackup/ummz/analyses/run_12_Aug20/5_counting_PE_noXY/arc_files/output.$JOB_ID.txt


# this line was different in Ian's script:    #$ -l h_rt=04:00:00

# (0) /home/home02/ummz/github_dirs/RNA-seq/scripts/run-5-counting.sh
# (1) running option [dups or nodups]
# (2) running mode parameter [SE or PE]
# (3) /path/to/input/data [bam files obtained from Picard] 
# (4) path to GTF annotation files [hg38.ncbiRefSeq.gtf]
# (5) /path/to/results 
# >> 
# (OUTPUT) /path/to/arc_files/output.$JOB_ID.txt

