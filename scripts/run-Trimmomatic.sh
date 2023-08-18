######################################################################################################
# this script runs trimming (on raw fastq files) for all files (pairs) at once
######################################################################################################

if [ $# != 4 ] ; then
    echo -e "ERROR: 4 arguments are required: \n
    (1) running mode (either 'SE' or 'PE'), \n
    (2) path to input data folder, \n
    (3) path to output folder and \n
    (4) SGE_TASK_ID argument for array jobs \n
    ... Exiting"
    exit 1
fi	

# define arguments
run_mode=$1         		# running mode (either 'SE' or 'PE')
data_dir=$2			# path to the folder with fastq files
out_dir=$3			# path the the folder where output will be placed

# export software (Trimmomatic)
export PATH=/nobackup/ummz/tools/bioinfo/Trimmomatic-0.39:$PATH

# get the read1 fastq.gz file and its pair
fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")

read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

coreFile=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run Trimmomatic
if [ $run_mode == 'SE' ]     # single-end [SE]
then
    echo "Running in single-end (SE) mode"
 
    java -jar /nobackup/ummz/tools/bioinfo/Trimmomatic-0.39/trimmomatic-0.39.jar \
    SE \
    -threads 4 \
    -phred33 \
    -trimlog ${coreFile}_trim_${run_mode}.log \
    -summary ${coreFile}_trim_${run_mode}_summary.txt \
    $read1 \
    $out_dir/${coreFile}_trim_${run_mode}_R1.fastq.gz \
    ILLUMINACLIP:/nobackup/ummz/tools/bioinfo/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
    SLIDINGWINDOW:4:15 \
    LEADING:20 \
    TRAILING:20 \
    HEADCROP:10 \
    MINLEN:36

    echo "Trimming completed"

elif [ $run_mode == 'PE' ]  # paired-end [PE]
then
    echo "Running in paired-end (PE) mode"
    
    java -jar /nobackup/ummz/tools/bioinfo/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE \
    -threads 4 \
    -phred33 -trimlog ${coreFile}_trim_${run_mode}.log -summary ${coreFile}_trim_${run_mode}_summary.txt \
    -validatePairs \
    $read1 $read2 \
    $out_dir/${coreFile}_trim_paired_R1.fastq.gz $out_dir/${coreFile}_trim_unpaired_R1.fastq.gz \
    $out_dir/${coreFile}_trim_paired_R2.fastq.gz $out_dir/${coreFile}_trim_unpaired_R2.fastq.gz \
    ILLUMINACLIP:/nobackup/ummz/tools/bioinfo/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
    SLIDINGWINDOW:4:15 \
    LEADING:20 \
    TRAILING:20 \
    HEADCROP:10 \
    MINLEN:36
    
    echo "Trimming completed"

else
    echo "ERROR: running mode argument incorrect! Exiting..."
    exit 1
fi
