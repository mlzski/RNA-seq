### this script run:wqs Salmon to obtaine reads from fastq files (trimmed) on transcript level
# NOTICE: prior to launching this command, make sure that the index has been already generated (just with a single command) 

# TODO: need to be tested for single-end data 

# this version is based on this tutorial: https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/salmon_kalliso_STAR_compare.md

# NOTE: make sure "Salmon" conda environment is activated !!!

if [ $# != 4 ] ; then
    echo -e "ERROR: 4 arguments are required: \
    (1) path to folder with input data, \
    (2) path to folder for output directory, \
    (3) path to index file \
    (4) SGE_TASK_ID argument for array jobs ...Exiting"
    exit 1
fi

# assign variables
data_dir=$1
out_dir=$2
index_dir=$3

# prepare params.txt file with information about the current run
me=$(basename "$0")
echo -ne "The same parameters were used for all samples." "\n" > params_${me}.txt
echo -ne "Script directory:" `pwd`"/"$me "\n" > params_${me}.txt

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
echo -ne "Executed on:" $current_date_time "\n" > params_${me}.txt

# get read1 fastq.gz file and its pair read2

fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")
read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

samp_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 25- | rev | sed -n -e "$SGE_TASK_ID p")

# run Salmon mapping job
salmon quant \
	-i $index_dir \
	-l A \
	-1 $read1 \
	-2 $read2 \
	-p 8 \
	--gcBias \
	--validateMappings \
	--writeUnmappedNames \
	--writeMappings \
	-o $out_dir/${samp_name}
 
# print the main command to params.txtSalmon_quantif_cc_c1_c2_transcript_lvl.o5169250.1
echo "salmon quant -i $index_dir -l A -1 $read1 -2 $read2 -p 8 --seqBias --gcBias --validateMappings --writeUnmappedNames --writeMappings -o $out_dir/${samp_name}" > params_${me}.txt

echo "Finished for " $samp_name 
