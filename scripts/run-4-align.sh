###################################################################################################### 
# this script runs STAR aligner to align the reads from fatq files to the genome
# there used to be 3 files: run-4-align_SE.sh, run-run-4-align_PE_paired.sh and run-4-align_PE_unpaired.sh, which were merged in this file
###################################################################################################### 

# change 'paired' to 'unpaired'

if [ $# != 4 ] ; then
    echo -e "ERROR: 4 arguments are required: \
    (1) path to data folder, \
    (2) path to output directory, \
    (3) path to index file and \
    (4) SGE_TASK_ID argument for array jobs \
    ...Exiting"
    exit 1
fi    

# export software (STAR)
export PATH=/home/home02/ummz/tools/STAR-2.7.3a/bin/Linux_x86_64_static:$PATH

# assign variables
run_mode=               # 'SE', 'PE-paired' or 'PE-unpaired'
data_dir=$1
out_dir=$2
index_dir=$3

# get the read1 fastq.gz file and its pair
if [ $run_mode == 'SE' ]
then
fastqFile=$(ls $data_dir/*_R1_single.fq | sed -n -e "$SGE_TASK_ID p")
read1=$fastqFile
bam_name=$(ls $data_dir/*_R1_single.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run STAR alignment in single-end mode [SE]
STAR \
--runMode alignReads \
--genomeDir $index_dir \
--runThreadN 8 \
--readFilesIn $read1 \
--outFileNamePrefix $out_dir/bam/${bam_name} \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--outSAMstrandField intronMotif

elif [ $run_mode == 'PE-paired' ]
then
fastqFile=$(ls $data_dir/*_R1_paired.fq | sed -n -e "$SGE_TASK_ID p")
read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')
bam_name=$(ls $data_dir/*_R1_paired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run STAR alignment in paired-end mode [PE-paired]
STAR \
--runMode alignReads \
--genomeDir $index_dir \
--runThreadN 8 \
--readFilesIn $read1 $read2 \
--outFileNamePrefix $out_dir/bam/${bam_name} \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--outSAMstrandField intronMotif

elif [ $run_mode == 'PE-unpaired' ]
fastqFile=$(ls $data_dir/*_R1_unpaired.fq | sed -n -e "$SGE_TASK_ID p")

read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

bam_name=$(ls $data_dir/*_R1_unpaired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run STAR alignment in paired-end mode [PE-unpaired]
STAR \
--runMode alignReads \
--genomeDir $index_dir \
--runThreadN 8 \
--readFilesIn $read1 $read2 \
--outFileNamePrefix $out_dir/bam/${bam_name} \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--outSAMstrandField intronMotif

else
echo ERROR... run_mode argument must be specified as: 'SE', 'PE-paired' or 'PE-unpaired'
fi
