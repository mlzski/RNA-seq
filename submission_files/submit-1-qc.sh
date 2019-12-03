# submission script for run_trimm_1.sh
# Michal Zulcinski 2019-11-13

# Run from the current directory with current environment
#$ -cwd -V

# ask for some time (hh:mm:ss max of 48:00:00)
#$ -l h_rt=8:00:00

# ask for some memory (by default, 1G, without a request)
#$ -l h_vmem=1.5G

# ask for 1 core
#$ -l np=1

# send emails when job starts and ends
#$ -m be
#$ -M @leeds.ac.uk

# now run the job
./run_trimm_1.sh

# probably need to redirect outpus to a new file

