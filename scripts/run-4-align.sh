### this script runs STAR aligner to align the reads from fatq files to the genome
# NOTICE: run 'index' separately | samtools is not included here

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
export PATH=/home/home02/ummz/tools/STAR-2.7.3a/bin/Linux_x86_64_static:$PATH

# assign variables
run_mode=$1               # 'SE' or 'PE'
data_dir=$2
out_dir=$3
index_dir=$4

if [ $run_mode == 'SE' ] ; then
    # get the read1 fastq.gz file
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

elif [ $run_mode == 'PE' ] ; then
    # get the read1 fastq.gz file and its pair read2
    fastqFile=$(ls $data_dir/*_R1_paired.fq | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    read2=$(echo $read1 | sed 's/R1/R2/g')
    
    bam_name=$(ls $data_dir/*_R1_paired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

    # run STAR alignment in paired-end mode [PE]
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
    echo "ERROR... run_mode argument must be specified as: 'SE' or 'PE'"
    exit 1
fi

# parameters meanings [for a mapping job]
# --runThreadN NumberOfThreads                      => defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node.
# --genomeDir /path/to/genomeDir                    => specifies path to the genome directory where genome indices where generated
# --readFilesIn /path/to/read1 [/path/to/read2]     => name(s) (with path) of the files containing the sequences to be mapped (e.g. RNA-seq FASTQ files); both read1 and read2 files have to be supplied for PE 
# --readFilesCommand UncompressionCommand           => uncompresses input files, where UncompressionCommand is the un-compression command that takes the file name as input parameter, and sends the uncompressed output to stdout. 
# --outSAMtype BAM SortedByCoordinate               => output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
# --outSAMattributes All                            => SAM attributes to be used 
# --outReadsUnmapped Fastx                          => will output unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads into separate file(s) Unmapped.out.mate1(2), formatted the same way as input read files (i.e. FASTQ or FASTA). 
