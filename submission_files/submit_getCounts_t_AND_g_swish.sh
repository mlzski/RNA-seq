# submission script to run get_counts_swish.R to obtain count matrix on gene- and transcript-level
# Michal Zulcinski 2021-04-15

#$ -cwd -V
#$ -l h_rt=01:00:00
#$ -l h_vmem=32G
#$ -m be
#$ -M ummz-arc-records@outlook.com

Rscript /home/home02/ummz/github_dirs/RNA-seq/scripts/submit_getCounts_t_AND_g_swish.sh "both-levels" /nobackup/ummz/analyses/run_17_Jun21/quants_all/transcript-level /nobackup/ummz/analyses/run_17_Jun21/swish/both-levels /nobackup/ummz/analyses/run_17_Jun21/swish/samples_all.txt
         
#############################################################
# Rscript run_DEG_swish_FINAL.R
# (1) INPUT directory
# (2) OUTPUT directory
# (3) annotation file [samples.txt]


