### this script runs second QC analysis (on trimmed fatstq files) for all files in parallel (using array task)

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: \
    (1) Path to data folder, \
    (2) path to output folder and \
    (3) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi	

# define arguments
data_dir=$1			# path to the folder with fastq files
out_dir=$2			# path the the folder where output will be placed

# export software (FastQC)
export PATH=/home/home02/ummz/tools/FastQC:$PATH

# check if output folders exist and create them if needed
if [ ! -d "$out_dir/report" ]; then
    echo -ne "\nCreating directories..."
    mkdir $out_dir/report $out_dir/temp
fi
 
# get the fastq.gz files names (no need to get read1 and read2 separately)
fastqFile=$(ls $data_dir/*.fq | sed -n -e "$SGE_TASK_ID p")

# run the fastqc command for each fastq.gz file
fastqc -o $out_dir/report --threads 4 --dir $out_dir/temp $fastqFile
