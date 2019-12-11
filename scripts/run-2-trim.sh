
################################################################################# 
# this script runs trimming (on raw fastq files) for all files (pairs) at once
#################################################################################

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: (1) Path to data folder, (2) path to output folder and (3) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi	

# define arguments
data_dir=$1			# path to the folder with fastq files
out_dir=$2			# path the the folder where output will be placed

# export software (TrimGalore and Python)
export PATH=/nobackup/ummz/tools/Trimmomatic-0.39:$PATH

# get the read1 fastq.gz file and its pair
fastqFile=$(ls $data_dir/*_R1_001.fastq.gz | sed -n -e "$SGE_TASK_ID p")

read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

coreFile=$(ls $data_dir/*_R1_001.fastq.gz | rev | cut -d '/' -f 1 | cut -c 16- | rev)

# run the trimmomatic command for each pair of fastq files
echo java -jar /nobackup/ummz/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 $read1 $read2 $out_dir/processed_fastq/${coreFile}R1_paired.fq $out_dir/processed_fastq/${coreFile}R1_unpaired.fq $out_dir/processed_fastq/${coreFile}R2_paired.fq $out_dir/processed_fastq/${coreFile}R2_unpaired.fq

