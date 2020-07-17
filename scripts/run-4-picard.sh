######################################################################################################
# this script runs Picard Tools to locate and tag duplicate reads in a BAM or SAM file
######################################################################################################

# need to create a new directory within /bam and place there all _Aligned.sortedByCoord.out.bam

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: \
	(1) path to data folder [BAM files], \
	(2) path to output directory, \
	(3) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi    

# assign variables
data_dir=$1   #                                                   | 11026_S12_Aligned.sortedByCoord.out.bam
out_dir=$2    # /nobackup/ummz/analyses/rerun_Ian/rerun_1/pic     | 11026_S12_Aligned.sortedByCoord.out.bam
#METRICS_FILE # /nobackup/ummz/analyses/rerun_Ian/rerun_1/pic     | 11026_S12_Aligned.sortedByCoord.out.bam.metrics.txt

#core_name=`echo $data_dir | rev | cut -d'/' -f 1 | rev`

bamFile=$(ls $data_dir/* | sed -n -e "$SGE_TASK_ID p")

core_name=$(ls $data_dir/* | rev | cut -d'/' -f 1 | rev | sed -n -e "$SGE_TASK_ID p")

java -jar /home/home02/ummz/tools/picard/build/libs/picard.jar MarkDuplicates \
	INPUT=$bamFile \
	OUTPUT=$out_dir/$core_name \
	METRICS_FILE=$out_dir/${core_name}.metrics.txt \
	TAGGING_POLICY=All \
	TMP_DIR=$out_dir/temp \
	VALIDATION_STRINGENCY=LENIENT \
	MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \
	SORTING_COLLECTION_SIZE_RATIO=0.25 \
	REMOVE_SEQUENCING_DUPLICATES=false \
	REMOVE_DUPLICATES=false \
	ASSUME_SORTED=false \
	DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
	PROGRAM_RECORD_ID=MarkDuplicates \
	PROGRAM_GROUP_NAME=MarkDuplicates \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
	VERBOSITY=INFO \
	QUIET=false \
	COMPRESSION_LEVEL=5 \
	MAX_RECORDS_IN_RAM=500000 \
	CREATE_INDEX=false \
	CREATE_MD5_FILE=false \
	GA4GH_CLIENT_SECRETS=client_secrets.json

