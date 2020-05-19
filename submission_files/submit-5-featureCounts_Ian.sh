# submission script for run-5-featureCounts_Ian.R
# Michal Zulcinski 2020-05-20

#$ -cwd -V
#$ -l h_rt=04:00:00
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -m be
#$ -M ummz-arc-records@outlook.com

Rscript /nobackup/ummz/analyses/rerun_Ian/rerun_3/scripts/run-5-featureCounts_Ian.R dups /nobackup/ummz/analyses/rerun_Ian/rerun_2/pic /nobackup/ummz/reference/Ian/hg38.gtf /nobackup/ummz/analyses/rerun_Ian/rerun_3/featCounts >> /nobackup/ummz/analyses/rerun_Ian/rerun_3/output.$JOB_ID.txt

# Rscript /nobackup/ummz/analyses/rerun_Ian/rerun_3/scripts/run-5-featureCounts_Ian.R nodups /nobackup/ummz/analyses/rerun_Ian/rerun_2/pic /nobackup/ummz/reference/Ian/hg38.gtf /nobackup/ummz/analyses/rerun_Ian/rerun_3/featCounts >> /nobackup/ummz/analyses/rerun_Ian/rerun_3/output.$JOB_ID.txt

