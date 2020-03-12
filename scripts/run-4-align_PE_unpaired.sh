###################################################################################################### 
# this script runs STAR aligner to align the reads from fatq files to the genome 
###################################################################################################### 

# change 'paired' to 'unpaired'

if [ $# != 4 ] ; then
    echo -e "ERROR: 4 arguments are required: (1) path to data folder, (2) path to output directory, (3) path to index file and (4) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi    

# export software (STAR)
export PATH=/home/home02/ummz/tools/STAR-2.7.3a/bin/Linux_x86_64_static:$PATH

# assign variables
data_dir=$1
out_dir=$2
index_dir=$3

# get the read1 fastq.gz file and its pair
fastqFile=$(ls $data_dir/*_R1_unpaired.fq | sed -n -e "$SGE_TASK_ID p")

read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

bam_name=$(ls $data_dir/*_R1_unpaired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run STAR alignment in paired-end mode [PE]
STAR --runMode alignReads --genomeDir $index_dir --runThreadN 8 --readFilesIn $read1 $read2 --outFileNamePrefix $out_dir/bam/${bam_name} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif


