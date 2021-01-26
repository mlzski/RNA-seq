### this script runs Salmon to obtaine reads from fastq files (trimmed) on transcript level
# NOTICE: prior to launching this command, make sure that the index has been already generated (just with a single command) 

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
export PATH=/home/home02/ummz/tools/salmon-latest_linux_x86_64/bin:$PATH

# assign variables
run_mode=$1               # 'SE' or 'PE'
data_dir=$2
out_dir=$3
index_dir=$4

# run Salmon mapping job
if [ $run_mode == 'SE' ] ; then
    # get the read1 fastq.gz file
    fastqFile=$(ls $data_dir/*_R1_single.fq | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    
    samp_name=$(ls $data_dir/*_R1_single.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

    # run Salmon in single-end mode [SE]
    #STAR \	[TODO: TO BE COMPLETED]
    #--runMode alignReads \
    #--genomeDir $index_dir \
    #--runThreadN 8 \
    #--readFilesIn $read1 \
    #--outFileNamePrefix $out_dir/bam/${bam_name} \
    #--outSAMtype BAM SortedByCoordinate \
    #--outSAMunmapped Within \
    #--outSAMattributes Standard 

elif [ $run_mode == 'PE' ] ; then
    # get the read1 fastq.gz file and its pair read2
    fastqFile=$(ls $data_dir/*_R1_paired.fq | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    read2=$(echo $read1 | sed 's/R1/R2/g')
    
    samp_name=$(ls $data_dir/*_R1_paired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

    # run Salmon in paired-end mode [PE]
    salmon quant \ 
	-i $index_dir \
	-l A \
	-1 $read1 \
	-2 $read2 \
	-p 8 \ 
	--validateMappings \ 
	-o $out_dir/${samp_name}_quant

else
    echo "ERROR... run_mode argument must be specified as: 'SE' or 'PE'"
    exit 1
fi

