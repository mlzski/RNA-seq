### this script runs SAMtools to index bam files from STAR aligner
# NOTE: to be run LOCALLY (no need for a submission file)

if [ $# != 3 ] ; then
    echo -e "ERROR: 2 argument are required: \
    (1) Path to data folder where _Aligned.sortedByCoord.out.bam file are stored, \
    (2) path to output directory and 
    (3) name for the .log file to be created \
    ...Exiting"
    exit 1
fi	

# USAGE EXAMPLE:
# bash run-4-samtools.sh /nobackup/ummz/analyses/rerun_FINAL/run_1/alignment_SE/bam /nobackup/ummz/analyses/rerun_FINAL/run_1 samtools_SE.log

# define arguments
data_dir=$1
out_dir=$2
log_file=$3

# export software (SAMtools)
export PATH=/home/home02/ummz/tools/samtools-1.10/bin:$PATH  

# go to data directory
cd $data_dir

# run SAMtools to index bam files
counter=0
for i in *_Aligned.sortedByCoord.out.bam 
do
    samtools index $i
    mv $i.bai $out_dir 
    echo Created $outdir/$i.bai  
    (( counter++ ))
done > $out_dir/$log_file

echo DONE. Indexed $counter files.
 
