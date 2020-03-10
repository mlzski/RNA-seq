
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
export PATH=/home/home02/ummz/tools/Trimmomatic-0.39:$PATH

# get the read1 fastq.gz file and its pair
fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")

read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

coreFile=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run the trimmomatic command [SE or PE] 
# NOTICE: adjust the parameters for each analysis

# single-end [SE]
java -jar /home/home02/ummz/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 -phred33 $read1 $out_dir/processed_fastq/${coreFile}R1_single.fq ILLUMINACLIP:/home/home02/ummz/tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:36

# paired-end [PE]
java -jar /home/home02/ummz/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 $read1 $read2 $out_dir/processed_fastq/${coreFile}R1_paired.fq $out_dir/processed_fastq/${coreFile}R1_unpaired.fq $out_dir/processed_fastq/${coreFile}R2_paired.fq $out_dir/processed_fastq/${coreFile}R2_unpaired.fq ILLUMINACLIP:/home/home02/ummz/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:36

