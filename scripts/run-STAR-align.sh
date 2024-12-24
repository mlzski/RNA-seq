### this script runs STAR aligner to align the reads from fastq files to the genome

# NOTE run 'index' separately | samtools is not included here

if [ $# != 6 ] ; then
    echo -e "ERROR: 6 arguments are required: \
    (1) running mode (either 'SE' or 'PE'), \
    (2) trimmed data? (yes/no), \
    (3) path to data folder, \
    (4) path to output directory, \
    (5) path to index file and \
    (6) SGE_TASK_ID argument for array jobs ...Exiting"
    exit 1
fi

# export software (STAR)
export PATH=/nobackup/ummz/tools/bioinfo/STAR-2.7.11a/bin/Linux_x86_64_static:$PATH

# assign variables
run_mode=$1               # 'SE' or 'PE'
trimmed=$2		  # 'yes' or 'no'
data_dir=$3
out_dir=$4
index_dir=$5

# run STAR alignment job
if [ $run_mode == 'SE' ] ; then

    # get .fastq.gz files for R1	

    fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    
    core_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -d '.' -f 3 | rev | sed 's/...$//' | sed -n -e "$SGE_TASK_ID p")
 
    if [ $trimmed == 'yes' ] ; then    
	bam_name=${core_name}
    elif [ $trimmed == 'no' ] ; then
	bam_name=${core_name}_no_trim
    fi
 
    mkdir -p $out_dir/bam/${core_name}_aligned

    # run STAR alignment in single-end mode [SE]
    STAR \
    --runMode alignReads \
    --genomeDir $index_dir \
    --runThreadN 8 \
    --readFilesCommand gunzip -c \
    --readFilesIn $read1 \
    --outFileNamePrefix $out_dir/bam/${core_name}_aligned/${bam_name}_ \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes All  
#    --outReadsUnmapped Fastx     # only for paired end
 

elif [ $run_mode == 'PE' ] ; then
    
    # get .fastq.gz files for R1 and their R2 pairs
    
    fastqFile=$(ls $data_dir/*_R1.fastq.gz | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    read2=$(echo $read1 | sed 's/R1/R2/g')
   
    core_name=$(ls $data_dir/*_R1.fastq.gz | rev | cut -d '/' -f 1 | cut -d '.' -f 3 | rev | sed 's/...$//' | sed -n -e "$SGE_TASK_ID p") 
    
    if [ $trimmed == 'yes' ] ; then
	bam_name=${core_name}
    elif [ $trimmed == 'no' ] ; then
	bam_name=${core_name}_no_trim
    fi

    mkdir -p $out_dir/bam/${core_name}_aligned

    # run STAR alignment in paired-end mode [PE]
    STAR \
    --runMode alignReads \
    --genomeDir $index_dir \
    --runThreadN 8 \
    --readFilesCommand gunzip -c \
    --readFilesIn $read1 $read2 \
    --outFileNamePrefix $out_dir/bam/${core_name}_aligned/${bam_name}_ \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes All \
    --outReadsUnmapped Fastx

else
    echo "ERROR... run_mode argument must be specified as: 'SE' or 'PE'"
    exit 1
fi

