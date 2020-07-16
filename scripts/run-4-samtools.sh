######################################################################################################
# this script runs SAMtools to index bam files from STAR aligner
######################################################################################################

if [ $# != 1 ] ; then
    echo -e "ERROR: 1 argument are required: \
    (1) Path to data folder where _Aligned.sortedByCoord.out.bam file are stored \
    ...Exiting"
    exit 1
fi	

# define arguments
$data_dir=$1

# export software (SAMtools)
export PATH=/home/home02/ummz/tools/samtools-1.10/bin:$PATH  

# go to data directory
cd $data_dir

# run SAMtools to index bam files
samtools index *_Aligned.sortedByCoord.out.bam
