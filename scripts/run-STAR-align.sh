### this script runs STAR aligner to align the reads from fastq files to the genome

# NOTE run 'index' separately | samtools is not included here

if [ $# != 5 ] ; then
    echo -e "ERROR: 5 arguments are required: \
    (1) running mode (either 'SE' or 'PE'), \
    (2) path to data folder, \
    (3) path to output directory, \
    (4) path to index file and \
    (5) SGE_TASK_ID argument for array jobs ...Exiting"
    exit 1
fi

# export software (STAR)
export PATH=/nobackup/ummz/tools/bioinfo/STAR-2.7.11a/bin/Linux_x86_64_static:$PATH

# assign variables
run_mode=$1               # 'SE' or 'PE'
data_dir=$2
out_dir=$3
index_dir=$4

# run STAR alignment job
if [ $run_mode == 'SE' ] ; then

    # get .fastq.gz files for R1	

    fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    
    bam_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")
    core_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 21- | rev | sed -n -e "$SGE_TASK_ID p")

    mkdir -p $out_dir/bam/${core_name}_SE_aligned

    # run STAR alignment in single-end mode [SE]
    STAR \
    --runMode alignReads \
    --genomeDir $index_dir \
    --runThreadN 8 \
    --readFilesCommand gunzip -c \
    --readFilesIn $read1 \
    --outFileNamePrefix $out_dir/bam/${core_name}_SE_aligned/${bam_name}_ \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes All  
#    --outReadsUnmapped Fastx     # only for paired end

elif [ $run_mode == 'PE' ] ; then
    
    # get .fastq.gz files for R1 and their R2 pairs
    
    fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    read2=$(echo $read1 | sed 's/R1/R2/g')
    
    bam_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")
    core_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -c 25- | rev | sed -n -e "$SGE_TASK_ID p")
    
    mkdir -p $out_dir/bam/${core_name}_PE_aligned

    # run STAR alignment in paired-end mode [PE]
    STAR \
    --runMode alignReads \
    --genomeDir $index_dir \
    --runThreadN 8 \
    --readFilesCommand gunzip -c \
    --readFilesIn $read1 $read2 \
    --outFileNamePrefix $out_dir/bam/${core_name}_PE_aligned/${bam_name}_ \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes All \
    --outReadsUnmapped Fastx

else
    echo "ERROR... run_mode argument must be specified as: 'SE' or 'PE'"
    exit 1
fi

